reads = open("./tmp.data/BAM/cis.txt",'r')
results = open("./tmp.data/BAM/result_inward.txt",'w')
results_outward = open("./tmp.data/BAM/result_outward.txt",'w')
results_samestrand = open("./tmp.data/BAM/result_samestrand.txt",'w')
inward=[([0] * 10000) for i in range(1)]
outward=[([0] * 10000) for i in range(1)]
samestrand=[([0] * 10000) for i in range(1)]
for read in reads:
	b  =  int(read.split(' ')[3]) - int(read.split(' ')[6])
	a  =  abs(b)
	if a < 100000:
		if int(read.split(' ')[1]) == int(read.split(' ')[4]):
			samestrand[0][a//10]+=1
		else:
			if read.split(' ')[1] == '1' and read.split(' ')[4] == '0':
				if b >0:inward[0][a//10]+=1
				else:   outward[0][a//10]+=1
			else:
				if b >0:outward[0][a//10]+=1
				else:   inward[0][a//10]+=1

f = str(inward).replace("], [","\t")
f = f.replace("[","")
f= f.replace("]","")
f= f.replace(","," ")
results.write(f+"\n")
results.close()


f = str(outward).replace("], [","\t")
f = f.replace("[","")
f= f.replace("]","")
f= f.replace(","," ")
results_outward.write(f+"\n")
results_outward.close()

f = str(samestrand).replace("], [","\t")
f = f.replace("[","")
f= f.replace("]","")
f= f.replace(","," ")
results_samestrand.write(f+"\n")
results_samestrand.close()
