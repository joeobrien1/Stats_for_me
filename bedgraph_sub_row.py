# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

#Read in input file Actual

testlist = dict()
with open("/Users/Kayleigh/Dropbox/Insight_Data_Engineer/input/actual.txt") as actualfile:
    for line in actualfile:
        parts = line.split(|)
        testlist.append = parts[0]
        

#Read in input file Prediction

#Read in input file Window
        
        
        
# Subset a bedgraph file

chrom='chr9'


new_bedgraph=[]
with open('/Users/Kayleigh/Desktop/Genome_mods/ENCFF118_H3K4me3.bedGraph','r') as bedgraph:
    for line in bedgraph:
        items=line.split()
        if chrom==items[0]:
            new_bedgraph.append(line)
bedgraph.close()


with open('/Users/Kayleigh/Desktop/Genome_mods/ENCFF208CCL_H2K27ac_sig2_chr9.bedGraph','w') as new_file:
    new_file.write('track type=bedGraph\n')
    for line in new_bedgraph:
        new_file.write(line)
new_file.close()

print('Completed')

