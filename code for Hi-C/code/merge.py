#!/usr/bin/python
# _*_ coding: UTF-8 _*_
import argparse
import sys
args=sys.argv
infilename=args[args.index('-i')+1]
outfilename=args[args.index('-o')+1]
samfile = open(infilename,'r')
result = open(outfilename,'w')
delimiter = ' '
a = samfile.readline()
for line in samfile:

    if line.split(' ')[0] == a.split(' ')[0]:
      # 将名字相同个的两个read写成一行：
        result.write(a.strip() +' '+ delimiter.join(line.split(' ')[1:4]))
    else: a = line
result.close()
