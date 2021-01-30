import argparse
import sys
args=sys.argv
infilename=args[args.index('-i')+1]
outfilename=args[args.index('-o')+1]
chrnum = args[args.index('-chr')+1]
resolution = args[args.index('-resolution')+1]
resolution = int(resolution)
reads = open(infilename,'r')
result = open(outfilename,'w')
chr = open("/lustre/user/liclab/lirf/Project/hic/bgi/5534N/code/chr.length.txt",'r')
a = chr.readlines()
num = int(a[int(chrnum)-1].strip().split('\t')[1])//resolution +1
data=[([0] * num) for i in range(num)]
for read in reads:
    read1 = read.split(' ')
    if int(read1[3])//resolution == int(read1[6])//resolution:
      data[int(read1[3])//resolution][int(read1[6])//resolution] += 1
    else:
      data[int(read1[3])//resolution][int(read1[6])//resolution] += 1
      data[int(read1[6])//resolution][int(read1[3])//resolution] += 1
f = str(data).replace("], [","\n")
f = f.replace("[","")
f= f.replace("]","")
f= f.replace(","," ")
result.write(f+"\n")
result.close()
