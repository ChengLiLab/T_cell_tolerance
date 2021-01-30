reads = open("./tmp.data/BAM/cis.txt",'r')
results = open("./tmp.data/BAM/sum.txt",'w')
inw = 0
out = 0
same = 0
for read in reads:
	b  =  int(read.split(' ')[3]) - int(read.split(' ')[6])
	a  =  abs(b)
	if int(read.split(' ')[1]) == int(read.split(' ')[4]):
			same+=1
	else:
			if read.split(' ')[1] == '1' and read.split(' ')[4] == '0':
				if b >0:inw+=1
				else:   out+=1
			else:
				if b >0:out+=1
				else:   inw+=1
results.write(str(inw)+"\n")
results.write(str(out)+"\n")
results.write(str(same)+"\n")
results.close()
