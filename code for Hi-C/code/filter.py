#!/usr/bin/python
# _*_ coding: UTF-8 _*_
### 此脚本用于将比对好的reads进行第一次过滤
import argparse
import sys
args=sys.argv
infilename=args[args.index('-i')+1]
outfilename=args[args.index('-o')+1]
samfile = open(infilename,'r')
result = open(outfilename,'w')
delimiter = ' '
delimiters = ''
for line in samfile:
    if 'Un_g' not in line:
        if 'hap'  not in line:
            if 'random' not in line:
                if 'chrM' not in line:
                    line = line.replace('chr','')
                    line = line.replace('X','23')
                    line = line.replace('Y','24')
                    if line.split(' ')[1] =='0':
                        result.write(line)
                    else:
                        result.write(line.split(' ')[0]+' 1 '+delimiter.join(line.split(' ')[2:4]))
result.close()
