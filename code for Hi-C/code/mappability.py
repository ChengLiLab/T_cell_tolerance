#!/usr/bin/python
# _*_ coding: UTF-8 _*_
# 此脚本用来计算每个bin 的mappability
# 方法：以酶切位点上下游各500bp长度为分析位置，以36bp长度，间隔9bp来构建sud0-reads: 55 * 2 = 110条reads （55 * 9 = 495bp)
from __future__ import division
import argparse
import sys
args=sys.argv
infilename=args[args.index('-i')+1]
raw=args[args.index('-raw')+1]
outfilename=args[args.index('-o')+1]
resolution=args[args.index('-resolution')+1]
chrnum=args[args.index('-chr')+1]
resolution = int(resolution)
result = open("./mappabilityForhindIII/"+outfilename,'w') # 此文件用于生成每个酶切位点的mappability

# 1. 产生sudo reads： 通过 ./code/creat_mappability_fa.py产生， /lustre/user/dengl/Project/HiC/BAM/result/chr1.fa
# 2. 用bowtie2 比对chr.fa ： -f chr1.fa -S chr1.sam ， -f 参数：比对fa文件，默认质量值相等，得到sam文件
# 3. sam文件过滤：awk -F '[ ;\t]' '$6 == "chr1"&& $5 == 0&& $4 == $7{print $1,$2,$3,$4,$5,$6,$7}' chr1.sam > chr1filter.txt
reads = open(infilename,'r') # 匹配正确的reads

# 4. 得到所有reads的名字信息： grep '>' chr1.fa > chr1.raw.txt
raw_locus = open(raw,'r') #所有的reads# 

# 5. 计算mappability 
chr = open("/lustre/user/liclab/lirf/Project/hic/bgi/5534N/code/chr.length.txt",'r')
a = chr.readlines()
num = int(a[int(chrnum)-1].strip().split('\t')[1])//resolution +1

data=[[0]*num] 

data1=[[0]*num]

for read in reads:
    data[0][int(float(read.split(' ')[1]))//resolution] += 1 # 计算每个bin 匹配正确的reads个数
    
for raw_locu in raw_locus:
    data1[0][int(float(raw_locu.strip().split(';')[3]))//resolution] += 1 # 计算每个bin中所有的reads个数
    
for i in range(num):
    if data1[0][i] != 0:
        result.write(str(data[0][i]/data1[0][i])+"\n") # 两者相除,即此bin的mappability
    else:result.write(str("0")+"\n")
    
result.close()

