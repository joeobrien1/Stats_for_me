#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 15:51:49 2020

@author: Kayleigh
"""

f=open('/Users/Kayleigh/Desktop/mm10_gene_info.txt','w+')

with open("/Users/Kayleigh/Desktop/Mus_musculus.GRCm38.100.chr.gtf",'r') as GTF_file:
    for line in GTF_file:
        if not line.startswith('#'):
            items=line.split('\t')
            if items[2] == 'gene':
                items2=items[8].split(';')
                gene_ID=items2[0][9:len(items2[0])-1]
                gene_name=items2[2][12:len(items2[2])-1]
                f.write('chr'+items[0]+'\t'+items[3]+'\t'+items[4]+'\t'+gene_ID+'\t'+gene_name+'\n')
                
            
    





