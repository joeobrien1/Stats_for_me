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
##### McNally Lab #####
#######################


#This is the basic version of the RNA-seq analysis script. It will take raw read files from paired-end RNAseq and generate gene counts, which can be inputted into other tools for differential expression analysis. 

#Important instructions
#1. Place all of your unzipped fastq files in an directory on quest. They must be named according to the following. One must end with ...R1_001.fastq and it's paired mate ends with ...R2_001.fastq. 
#2. Choose a reference. hg19 or mm10 are supported
#3. Run this python script in quest with the options. Use -h for options
#4. Submit scripts to quest for running
#5. Wait and relax

#Supported options
#- Can use STAR or TopHat for alignment
#- supportes hg19 and mm10 references
#- Multiple stranded options. The mRNA Truseq kit that NUseq uses requires reverse strand

##Suggested Steps:
# - FastQC, Trim, and Align with STAR, HT-seq Count

##SETUP##
import os
from optparse import OptionParser

#Command line Arguments#
parser=OptionParser()
parser.add_option("-w","--working_directory", dest="working_directory", default=" ", help="Put the full working directory path here (where the fastq files are located) Default is current directory.")
parser.add_option("-t","--type", dest="analysis_type", default="counts", help="What type of analysis do you want? Default= counts")
parser.add_option("-a","--aligner", dest="aligner", default="STAR", help="STAR, tophat, none, or both. Default=STAR")
parser.add_option("-q","--fastQC", dest="fastQC", default='1', help="Run FastQC? Default= 1")
parser.add_option("-x","--trim", dest="trim", default='1', help="Do you want to trim the reads? Default= 1")
parser.add_option("-r","--reference", dest="reference", default=' ', help="Specify reference genome to use. hg19 and mm10 are supported. Required parameter.")
parser.add_option("-c","--strand", dest="strand", default='reverse', help="Is your RNA-seq data stranded? Options: yes, no, reverse. Default=reverse")
parser.add_option("-u","--user", dest="user", default =" ", help= "Denote a user ID. Net_ID work well. Required parameter.")
parser.add_option("-g","--group", dest="allocation", default ="b1042", help= "Denote an allocation number. Default is the genomics cluster b1042")
parser.add_option("-p","--paired",dest="paired", default="1", help= "Is your data paired end or single end sequncing? Options are 1 or 0. Default 1")

(options,args) = parser.parse_args()

# Tools #may need to move these from a permissions standpoint
Trimmer= '/projects/genomicsshare/rna-seq_reference/McNally_TOOLS/trimmomatic-0.36.jar'

#Reference Files (hg19) 
#Reference files are stored on quest at these shared directories
if options.reference == 'hg19':
    Ensembl_GTF_Tophat='/projects/genomicsshare/rna-seq_reference/McNally_Reference/Ensembl/known'
    Ensembl_GTF='/projects/genomicsshare/rna-seq_reference/McNally_Reference/Ensembl/Homo_sapiens.GRCh37.87.gtf'
    Ensembl_bowtie_index='/projects/genomicsshare/rna-seq_reference/McNally_Reference/Ensembl/Homo_sapiens.GRCh37'
    Ensembl_Genome_Dir='/projects/genomicsshare/rna-seq_reference/McNally_Reference/Ensembl/GenomeDirectory'

#Reference Files (mm10)
#Reference files are stored on quest at these directories
if options.reference=='mm10':
    Ensembl_GTF_Tophat='/projects/p20742/anno/tophat_tx/mm10.Ens_78.cuff'
    Ensembl_GTF='/projects/p20742/anno/Ens/mm10.Ens_78/mm10.Ens_78.cuff.gtf '
    Ensembl_bowtie_index='/projects/p20742/anno/tophat_tx/mm10.Ens_78.cuff'
    Ensembl_Genome_Dir='/projects/p20742/anno/STAR_indexes/mm10/'


##IGNORE
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

if options.reference == ' ':
    print('No reference indicated, please select hg19 or mm10 with the -r tag')
    quit()
if options.user == ' ':
    print('No user indicated, please add your netID with the -u tag')
    quit()
    
    
InputDir = options.working_directory
listdir = os.listdir(InputDir)

#create folders if they dont already exist
RootFolder = InputDir+'/RNA_Seq_Analysis_'+options.user+'/'
ScriptFolder = InputDir+'/RNA_Seq_Analysis_'+options.user+'/Scripts'
logFolder= InputDir+'/RNA_Seq_Analysis_'+options.user+'/logs'


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
    - working_directory =  '''+options.working_directory+'''
    - FastQC = '''+options.fastQC+'''
    - trim = '''+options.trim+'''
    - strand = '''+options.strand+'''
    - user = '''+options.user+'''
    - allocation = '''+options.allocation+'''
    - paried = '''+options.paired+'''
-------------------------------------    
''')    
    
for file in listdir:
    if file.endswith("R1_001.fastq"):
        #Makes Shell script with file name
        script = open(ScriptFolder+'/'+file[0:len(file)-12]+'_Master_Run.sh', "w+")
        #Print MSUB Header
        script.write('''#!/bin/bash
#SBATCH -A '''+options.allocation+'''
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
            R1Path = Trimmed_Dir+'/'+file[0:len(file)-12]+'_F_paired.fastq'
            R2Path = Trimmed_Dir+'/'+file[0:len(file)-12]+'_R_paired.fastq'
        
        #Align with Tophat
        if options.paired=='1':
            Tophat_Dir= FileFolder+'/Tophat'
            if options.aligner == 'tophat' or options.aligner == 'both':                                                                                        
                print >> script, 'mkdir '+Tophat_Dir
                print >> script, 'tophat --no-novel-juncs --transcriptome-index '+Ensembl_GTF_Tophat+' --num-threads 6 --max-multihits 5 -o '+Tophat_Dir+' '+Ensembl_bowtie_index+' '+R1Path+' '+R2Path
                print >> script, 'samtools sort -o '+Tophat_Dir+'/accepted_hits_sorted.bam '+Tophat_Dir+'/accepted_hits.bam'
                print >> script, 'samtools index -b '+Tophat_Dir+'/accepted_hits_sorted.bam' 
                
            #Align with STAR
            STAR_Dir= FileFolder+'/STAR_Mapping'
            if options.aligner == 'STAR' or options.aligner == 'both':
                print >> script, 'mkdir '+STAR_Dir                                                                                     
                print >> script, 'cd '+STAR_Dir
                print >> script, 'STAR --genomeDir '+Ensembl_Genome_Dir+' --readFilesIn '+R1Path+' '+R2Path+' --runThreadN 16 --twopassMode Basic --twopass1readsN -1 --outSAMtype BAM SortedByCoordinate --outSJfilterReads Unique '
                print >> script, 'cd '+FileFolder
        
        if options.paired=='0':
            Tophat_Dir= FileFolder+'/Tophat'
            if options.aligner == 'tophat' or options.aligner == 'both':                                                                                        
                print >> script, 'mkdir '+Tophat_Dir
                print >> script, 'tophat --no-novel-juncs --transcriptome-index '+Ensembl_GTF_Tophat+' --num-threads 6 --max-multihits 5 -o '+Tophat_Dir+' '+Ensembl_bowtie_index+' '+R1Path
                print >> script, 'samtools sort -o '+Tophat_Dir+'/accepted_hits_sorted.bam '+Tophat_Dir+'/accepted_hits.bam'
                print >> script, 'samtools index -b '+Tophat_Dir+'/accepted_hits_sorted.bam' 
                
            #Align with STAR
            STAR_Dir= FileFolder+'/STAR_Mapping'
            if options.aligner == 'STAR' or options.aligner == 'both':
                print >> script, 'mkdir '+STAR_Dir                                                                                     
                print >> script, 'cd '+STAR_Dir
                print >> script, 'STAR --genomeDir '+Ensembl_Genome_Dir+' --readFilesIn '+R1Path+' --runThreadN 16 --twopassMode Basic --twopass1readsN -1 --outSAMtype BAM SortedByCoordinate --outSJfilterReads Unique '
                print >> script, 'cd '+FileFolder
     
        #HT-Seq Count
        HTSeq= FileFolder+'/HT-Seq_Counts'
        if options.analysis_type == 'counts' or options.analysis_type == 'both':
            print >> script, 'mkdir '+HTSeq
            #for differential expression
            if options.aligner == 'STAR' or options.aligner == 'both':
                print >> script, 'htseq-count -q --order=pos -s '+options.strand+' -t exon -m intersection-nonempty --format=bam --idattr=gene_id --additional-attr=gene_name '+STAR_Dir+'/Aligned.sortedByCoord.out.bam '+Ensembl_GTF+' > '+HTSeq+'/'+file[0:len(file)-12]+'STAR_counts_diffexpr.txt'
            if options.aligner == 'tophat' or options.aligner == 'both':
                print >> script, 'htseq-count -q --order=pos -s '+options.strand+' -t exon -m intersection-nonempty --format=bam --idattr=gene_id --additional-attr=gene_name '+Tophat_Dir+'/accepted_hits_sorted.bam '+Ensembl_GTF+' > '+HTSeq+'/'+file[0:len(file)-12]+'Tophat_counts_diffexpr.txt'            

                                                                                                                                                                                                                                                        
        print >> script, 'echo \'Analysis completed for '+file[0:len(file)-12]+'\''
        script.close()        
        print ('Scripts created successfully for '+file[0:len(file)-12]+'')                                                                                           

print 'To run all scripts: find -name *__Master_Run.sh -exec sbatch {} \;'



