#!/usr/bin/python
# _*_ coding: UTF-8 _*_
import argparse
import sys
args=sys.argv
infilename=args[args.index('-i')+1]
outfilename=args[args.index('-o')+1]
genome = open('/lustre/user/liclab/lirf/Project/hic/hg19.filter2.txt','r')
reads = open(infilename,'r')
result = open(outfilename,'w')
data=[]
for chrom in genome:
        a = chrom.strip().upper()
        data.append(a)
for read in reads:
        n = int(read.split(' ')[3])
        if read.split(' ')[1] == '0':
            f = data[int(read.split(' ')[2])-1][n:(n+536)] # all the alignment start location locates at the 5' of the positive chain
            if f.find('GATC') != -1:
                    result.write(read)
        else:
            f = data[int(read.split(' ')[2])-1][(n-500):n]
            if f.find('GATC') != -1:
                result.write(read)
result.close()

