#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 26 14:06:16 2020

@author: Kayleigh
"""

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 11 17:01:25 2017

@author: amg171
"""
#######################
###RNA Seq Analysis####
#######################
#### Anthony Gacita ###
#### McNally Lab ######
#######################

#Three modes here
#1. Tophat, RNA-SeQC, compare with GTEX
#2. STAR, HT-Seq counts 
#3. Alternative Splicing Analysis with MISO

##Steps:
# 1. FastQC, Trim, and Align with STAR
# 2. Process BAM File for GATK
# 3. Coverage Plots- per gene 
# 4. Get counts with RNA-SeQC
# 5. Compare RNA-SeQC counts with GTEX
# 6. Get counts with HT-Seq- may need to test stranded parameter




##SETUP##
import os
from optparse import OptionParser

#Command line Arguments#
parser=OptionParser()
parser.add_option("-w","--working_directory", dest="working_directory", default=" ", help="Put the full working directory path here (where the fastq files are located) Default is current directory. No last /")
parser.add_option("-t","--type", dest="analysis_type", default="both", help="What type of analysis do you want? Set to gtex, counts, or both. Default= both")
parser.add_option("-a","--aligner", dest="aligner", default="both", help="STAR, tophat, none, or both. Default=Both")
parser.add_option("-o","--output", dest="output", default="./RNA_Seq_Analysis_AG/Scripts/", help="Where do you want to write the scripts? Default= Working_Directory/RNA_Seq_Analysis_AG/Scripts/")
parser.add_option("-q","--fastQC", dest="fastQC", default='1', help="Run FastQC? Default= 1")
parser.add_option("-x","--trim", dest="trim", default='1', help="Do you want to trim the reads? Default= 1")
parser.add_option("-s","--splicing", dest="splicing", default='1', help="Do you want to do alternative splicing anlysis? Default= 1")
parser.add_option("-r","--reference", dest="reference", default='-', help="Specify reference genome to use. hg19 and mm10 are supported")
parser.add_option("-c","--strand", dest="strand", default='reverse', help="Is your RNA-seq data stranded? Options: yes, no, reverse. Default=reverse")

(options,args) = parser.parse_args()

# Tools
picard = '/home/amg171/McNally/TOOLS/picard.jar'
GATK = '/home/amg171/McNally/TOOLS/GenomeAnalysisTK.jar'
Trimmer= '/home/amg171/McNally/TOOLS/trimmomatic-0.36.jar'
RNA_SeQC='/home/amg171/McNally/TOOLS/RNA-SeQC_v1.1.8.jar'

#Alt Splicing Settings
Read_length_splicing= '100'
MISO_event_choices=['SE','RI','MXE','A3SS','A5SS','isoforms']
pack_miso= 0

#Reference Files (hg19)
if options.reference == 'hg19':
    Ensembl_GTF_Tophat='/home/amg171/McNally/ensembl_GTF_index/known'
    Ensembl_GTF='/projects/genomicsshare/rna-seq_reference/McNally_Reference/Ensembl/Homo_sapiens.GRCh37.87.gtf'
    Ensembl_bowtie_index='/home/amg171/McNally/ensembl_bowtie2_index/Homo_sapiens.GRCh37'
    Ensembl_Genome_Dir='/projects/genomicsshare/rna-seq_reference/McNally_Reference/Ensembl/GenomeDirectory'
    Ensembl_Reference='/projects/genomicsshare/rna-seq_reference/McNally_Reference/Ensembl/Homo_sapiens.GRCh37.dna.primary_assembly.fasta'
    
    
    Gtex_GTF_Tophat='/projects/genomicsshare/rna-seq_reference/McNally_Reference/Gtex_References/transcriptome_data/known'
    Gtex_GTF='/projects/genomicsshare/rna-seq_reference/McNally_Reference/Gtex_References/gencode.v19.genes.v6p_model.patched_contigs.gtf'
    Gtex_Bowtie_Index='/projects/genomicsshare/rna-seq_reference/McNally_Reference/Gtex_References/Homo_sapiens_assembly19'
    Gtex_Genome_Dir='/projects/genomicsshare/rna-seq_reference/McNally_Reference/Gtex_References/GenomeDirectory'
    Gtex_Reference='/projects/genomicsshare/rna-seq_reference/McNally_Reference/Gtex_References/Homo_sapiens_assembly19.fa'
    Gtex_Reference_Expression='/projects/genomicsshare/rna-seq_reference/McNally_Reference/Gtex_References/GTEX_median_LV.gct'
    MISO_const_exons= '/projects/genomicsshare/rna-seq_reference/McNally_Reference/Ensembl/MISO/ensembl_miso_files/Homo_sapiens.GRCh37.87.min_1000.const_exons.gff'
    MISO_index='/projects/genomicsshare/rna-seq_reference/McNally_Reference/Ensembl/MISO/ensembl_miso_files/indexed' 


#Reference Files (mm10)
if options.reference=='mm10':
    Ensembl_GTF_Tophat='/projects/p20742/anno/tophat_tx/mm10.Ens_78.cuff'
    Ensembl_GTF='/projects/p20742/anno/Ens/mm10.Ens_78/mm10.Ens_78.cuff.gtf '
    Ensembl_bowtie_index='/projects/p20742/anno/tophat_tx/mm10.Ens_78.cuff'
    Ensembl_Genome_Dir='/projects/p20742/anno/STAR_indexes/mm10/'
    Ensembl_Reference=''



#MISO_gff= '/projects/genomicsshare/rna-seq_reference/McNally_Reference/Ensembl/MISO/hg19/SE.hg19.no_chr.gff3'
#MISO_const_exons= '/projects/genomicsshare/rna-seq_reference/McNally_Reference/Ensembl/MISO/ensembl_miso_files/Homo_sapiens.GRCh37.87.min_1000.const_exons.gff'
#MISO_index_exons= '/projects/b1042/McNallyLab/Anthony/Reference/MISO/indexed_SE_exons'
#MISO_index_isoforms='/projects/b1042/McNallyLab/Anthony/Reference/MISO/indexed_SE_isoforms'
#MISO_index='/projects/genomicsshare/rna-seq_reference/McNally_Reference/Ensembl/MISO/ensembl_miso_files/indexed' 
#MISO_index='/projects/b1042/McNallyLab/Anthony/Reference/MISO'


##SETUP END##


######### Script Start #########

if options.working_directory ==' ':
    options.working_directory=os.getcwd()

InputDir = options.working_directory
listdir = os.listdir(InputDir)

#create folders if they dont already exist
RootFolder = InputDir+'/RNA_Seq_Analysis_AG'
ScriptFolder = InputDir+'/RNA_Seq_Analysis_AG/Scripts'
logFolder= InputDir+'/RNA_Seq_Analysis_AG/logs'


if not os.path.exists(RootFolder):
    os.makedirs(RootFolder)
    
if not os.path.exists(ScriptFolder):
    os.makedirs(ScriptFolder)

if not os.path.exists(logFolder):
    os.makedirs(logFolder)        
    


print ('''
-------------------------------------
---RNA-Seq Analysis Python Script----
------Created by Anthony Gacita------
-------------------------------------

-Options:
    - analysis_type = '''+options.analysis_type+'''
    - aligner = '''+options.aligner+'''
    - script_directory = '''+options.output+'''
    - working_directory =  '''+options.working_directory+'''
    - FastQC = '''+options.fastQC+'''
    - trim = '''+options.trim+'''
    - splicing= '''+options.splicing+'''
    - strand = '''+options.strand+'''
-------------------------------------    
''')    
    
    
for file in listdir:
    if file.endswith("R1_001.fastq"):
        #Makes Shell script with file name
        script = open(InputDir+'/RNA_Seq_Analysis_AG/Scripts/'+file[0:len(file)-12]+'_Master_Run.sh', "w+")
        #Print MSUB Header
        script.write('''#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 48:00:00
#SBATCH --mail-type=BEGIN,END,NONE,FAIL,REQUEUE
#SBATCH --output='''+logFolder+'''/'''+file[0:len(file)-12]+'''_log
#SBATCH -N 1
#SBATCH --ntasks-per-node=8
#SBATCH --mem 50000
module load bowtie2/2.2.6
module load tophat/2.1.0
module load samtools
module load bedtools
module load boost
module load bowtie
module load gcc/4.8.3
module load java/jdk1.7.0_04 
module load R/3.3.1
module load STAR/2.5.2 
module load cufflinks
module load python
module load fastqc


''')      


        #Set Up File-Specific Folders
        FileFolder = RootFolder+'/'+file[0:len(file)-12]
        print >> script, 'mkdir '+FileFolder
        R1Path = InputDir+'/'+file
        R2Path = InputDir+'/'+file[0:len(file)-11]+'2_001.fastq'
        
        #FastQC
        QC_Dir= FileFolder+'/FastQC'
        if options.fastQC== '1':
            print >> script, 'mkdir '+QC_Dir
            print >> script, 'fastqc -t 2 -o '+QC_Dir+' '+R1Path+' '+R2Path
        
        #Trim the Reads
        Trimmed_Dir= FileFolder+'/TrimmedReads'
        if options.trim== '1':
            print >> script, 'mkdir '+Trimmed_Dir
            print >> script, 'java -jar '+Trimmer+' PE '+R1Path+' '+R2Path+' '+Trimmed_Dir+'/'+file[0:len(file)-12]+'_F_paired.fastq '+Trimmed_Dir+'/'+file[0:len(file)-12]+'_F_unpaired.fastq '+Trimmed_Dir+'/'+file[0:len(file)-12]+'_R_paired.fastq '+Trimmed_Dir+'/'+file[0:len(file)-12]+'_R_unpaired.fastq ILLUMINACLIP:/projects/genomicsshare/rna-seq_reference/McNally_Reference/adapters.fa:2:30:10 TRAILING:30 MINLEN:50'

        #Align with Tophat
        Tophat_Dir= FileFolder+'/Tophat'
        if options.aligner == 'tophat' or options.aligner == 'both':                                                                                        
            print >> script, 'mkdir '+Tophat_Dir
            print >> script, 'tophat --no-novel-juncs --transcriptome-index '+Ensembl_GTF_Tophat+' --num-threads 6 --max-multihits 5 -o '+Tophat_Dir+' '+Ensembl_bowtie_index+' '+Trimmed_Dir+'/'+file[0:len(file)-12]+'_F_paired.fastq '+Trimmed_Dir+'/'+file[0:len(file)-12]+'_R_paired.fastq' 
            print >> script, 'samtools sort -o '+Tophat_Dir+'/accepted_hits_sorted.bam '+Tophat_Dir+'/accepted_hits.bam'
            print >> script, 'samtools index -b '+Tophat_Dir+'/accepted_hits_sorted.bam' 
            
        #Align with STAR
        STAR_Dir= FileFolder+'/STAR_Mapping'
        if options.aligner == 'STAR' or options.aligner == 'both':
            print >> script, 'mkdir '+STAR_Dir                                                                                     
            print >> script, 'cd '+STAR_Dir
            print >> script, 'STAR --genomeDir '+Ensembl_Genome_Dir+' --readFilesIn '+Trimmed_Dir+'/'+file[0:len(file)-12]+'_F_paired.fastq '+Trimmed_Dir+'/'+file[0:len(file)-12]+'_R_paired.fastq --runThreadN 16 --twopassMode Basic --twopass1readsN -1 --outSAMtype BAM SortedByCoordinate --outSJfilterReads Unique '
            print >> script, 'cd '+FileFolder
                                                                                                    
        #GTEX Analysis
        Gtex_Analysis= FileFolder+'/GTEX_analysis'
        if options.analysis_type == 'gtex' or options.analysis_type == 'both':
            print >> script, 'mkdir '+Gtex_Analysis
            print >> script, 'mkdir '+Gtex_Analysis+'/Tophat_Alignment'
            print >> script, 'tophat --no-novel-juncs --transcriptome-index '+Gtex_GTF_Tophat+' --num-threads 6 --max-multihits 5 -o '+Gtex_Analysis+'/Tophat_Alignment/ '+Gtex_Bowtie_Index+' '+Trimmed_Dir+'/'+file[0:len(file)-12]+'_F_paired.fastq '+Trimmed_Dir+'/'+file[0:len(file)-12]+'_R_paired.fastq'                                                                     
            print >> script, 'java -jar '+picard+' AddOrReplaceReadGroups I='+Gtex_Analysis+'/Tophat_Alignment/accepted_hits.bam  O='+Gtex_Analysis+'/'+file[0:len(file)-12]+'AddedRG.bam  RGID=id1 RGLB=library RGPL=Illumina RGPU=machine1 RGSM=sample1' 
            print >> script, 'java -jar '+picard+' ReorderSam I='+Gtex_Analysis+'/'+file[0:len(file)-12]+'AddedRG.bam O='+Gtex_Analysis+'/'+file[0:len(file)-12]+'reordered.bam R='+Gtex_Reference
            print >> script, 'java -jar '+picard+' MarkDuplicates I='+Gtex_Analysis +'/'+file[0:len(file)-12]+'reordered.bam O='+Gtex_Analysis+'/'+file[0:len(file)-12]+'_RG_dups.bam  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M='+Gtex_Analysis+'/output.metrics'
            print >> script, 'rm '+Gtex_Analysis+'/'+file[0:len(file)-12]+'AddedRG.bam '+Gtex_Analysis+'/'+file[0:len(file)-12]+'reordered.bam' 
            #RNA-SeQC
            RNASeQC= Gtex_Analysis+'/RNA-SeQC'
            print >> script, 'mkdir '+RNASeQC
            print >> script, 'module unload java/jdk1.8.0_25' 
            print >> script, 'module load java/jdk1.7.0_04 '
            print >> script, 'java -jar '+RNA_SeQC+' -r '+Gtex_Reference+' -o '+RNASeQC+'/ -t '+Gtex_GTF+' -corr '+Gtex_Reference_Expression+' -s "'+file[0:len(file)-12]+'|'+Gtex_Analysis+'/'+file[0:len(file)-12]+'_RG_dups.bam|Notes" -strictMode'   
        
        #HT-Seq Count
        HTSeq= FileFolder+'/HT-Seq_Counts'
        if options.analysis_type == 'counts' or options.analysis_type == 'both':
            print >> script, 'mkdir '+HTSeq
            #for differential expression
            if options.aligner == 'STAR' or options.aligner == 'both':
                print >> script, 'htseq-count -q --order=pos -s '+options.strand+' -t exon -m intersection-nonempty --format=bam --idattr=gene_id --additional-attr=gene_name '+STAR_Dir+'/Aligned.sortedByCoord.out.bam '+Ensembl_GTF+' > '+HTSeq+'/'+file[0:len(file)-12]+'STAR_counts_diffexpr.txt'
            if options.aligner == 'tophat' or options.aligner == 'both':
                print >> script, 'htseq-count -q --order=pos -s '+options.strand+' -t exon -m intersection-nonempty --format=bam --idattr=gene_id --additional-attr=gene_name '+Tophat_Dir+'/accepted_hits_sorted.bam '+Ensembl_GTF+' > '+HTSeq+'/'+file[0:len(file)-12]+'Tophat_counts_diffexpr.txt'            
            #for alternative splicing
            #print >> script, 'htseq-count -q --order=pos -s reverse -t gene -m union --nonunique=all --format=bam --idattr=gene_id --additional-attr=gene_name -a 255 '+STAR_Dir+'/Aligned.sortedByCoord.out.bam '+Ensembl_GTF+' > '+HTSeq+'/'+file[0:len(file)-12]+'STAR_counts_splicing.txt'
        
        #EDGEr Analysis
        #EdgeR= FileFolder+'/EdgeR_Diff_Expression'
        #print >> script, 'mkdir '+EdgeR


        #MISO Splicing Analysis 
        #compute insert length, run miso, summarize                                                                                                                                                                                                                                            
        Alternative_Splicing= FileFolder+'/Alt_splicing'
        if options.splicing == '1':
            print >> script, 'mkdir '+Alternative_Splicing
            print >> script, 'java -jar '+Trimmer+' PE '+R1Path+' '+R2Path+' '+Alternative_Splicing+'/'+file[0:len(file)-12]+'_F_paired_samelength.fastq '+Alternative_Splicing+'/'+file[0:len(file)-12]+'_F_unpaired_samelength.fastq '+Alternative_Splicing+'/'+file[0:len(file)-12]+'_R_paired_samelength.fastq '+Alternative_Splicing+'/'+file[0:len(file)-12]+'_R_unpaired_samelength.fastq ILLUMINACLIP:/projects/genomicsshare/rna-seq_reference/McNally_Reference/adapters.fa:2:30:10 CROP:'+Read_length_splicing+' MINLEN:'+Read_length_splicing
            print >> script, 'mkdir '+Alternative_Splicing+'/Alignment/'
            print >> script, 'cd '+Alternative_Splicing+'/Alignment/'
            print >> script, 'STAR --genomeDir '+Ensembl_Genome_Dir+' --readFilesIn '+Alternative_Splicing+'/'+file[0:len(file)-12]+'_F_paired_samelength.fastq '+Alternative_Splicing+'/'+file[0:len(file)-12]+'_R_paired_samelength.fastq --runThreadN 16 --twopassMode Basic --twopass1readsN -1 --outSAMtype BAM SortedByCoordinate --outSJfilterReads Unique '
            print >> script, 'cd '+Alternative_Splicing
            print >> script, 'samtools index -b '+Alternative_Splicing+'/Alignment/Aligned.sortedByCoord.out.bam'
            print >> script, 'pe_utils --compute-insert-len '+Alternative_Splicing+'/Alignment/Aligned.sortedByCoord.out.bam '+MISO_const_exons+' --output-dir '+Alternative_Splicing+'/insert-dist/ --no-bam-filter'                                                                                                                                                                                                                                  
            print >> script, '''LINE1=$(awk 'FNR==1' '''+Alternative_Splicing+'''/insert-dist/Aligned.sortedByCoord.out.bam.insert_len) \nmean=${LINE1:6:5} \necho "Mean insert = $mean"'''
            print >> script, '''LINE1=$(awk 'FNR==1' '''+Alternative_Splicing+'''/insert-dist/Aligned.sortedByCoord.out.bam.insert_len) \nsdev=${LINE1:17:4} \necho "Stdev = $sdev"'''
            for choice in MISO_event_choices:
                print >> script, 'miso --run '+MISO_index+'/'+choice+' '+Alternative_Splicing+'/Alignment/Aligned.sortedByCoord.out.bam --output-dir '+Alternative_Splicing+'/'+file[0:len(file)-12]+'_MISO_out_'+choice+'/ --read-len '+Read_length_splicing+' --paired-end $mean $sdev'
                print >> script, 'summarize_miso --summarize-samples '+Alternative_Splicing+'/'+file[0:len(file)-12]+'_MISO_out_'+choice+'/ '+Alternative_Splicing+'/'+file[0:len(file)-12]+'/summary_'+choice                                                                                                                                                                                       
                if pack_miso == 1:
                    print >> script, 'miso_pack --pack '+Alternative_Splicing+'/'+file[0:len(file)-12]+'_MISO_out_'+choice+'/'
                                                                                                                                                                                                                                                        
        print >> script, 'echo \'Analysis completed for '+file[0:len(file)-12]+'\''
        script.close()        
        print ('Scripts created successfully for '+file[0:len(file)-12]+'')                                                                                           

print 'To run all scripts: find -name *__Master_Run.sh -exec sbatch {} \;'