a= read.table("./tmp.data/BAM/result_inward.txt",header =F)
b= read.table("./tmp.data/BAM/result_outward.txt",header =F)
c= read.table("./tmp.data/BAM/result_samestrand.txt",header =F)
sum.num =  read.table("./tmp.data/BAM/sum.txt",header =F)
a = as.matrix(a)
b = as.matrix(b)
c = as.matrix(c)
sum.a = sum.num[1,1]
sum.b = sum.num[2,1]
sum.c = sum.num[3,1]
sum = sum.a + sum.b + sum.c
par(cex=2)
pdf("./Graph/Cutoff distance between read pairs.pdf",5,5)
plot(0,xlim = c(1,2000),ylim = c(0,60),col="white",ylab = "Percentage of total reads (%)",
     xlab="Cutoff distance between read pairs ")
f = array(dim=c(196,4))
for(i in 5:200){
f[i-4,1] = i*10
f[i-4,2] = ((sum.a - sum(a[1:i]))/sum)*100
f[i-4,3] = ((sum.b - sum(b[1:i]))/sum)*100
f[i-4,4] = ((sum.c - sum(c[1:i]))/sum)*100
}
lines(f[,1],f[,2],type="l",col =1,lwd=6)
lines(f[,1],f[,3],type="l",col =2,lwd=6)
lines(f[,1],f[,4],type="l",col =3,lwd=6)
legend("topright", legend=c("inward","outward","samestrand"), col=1:3, lwd=2,border ="white",cex = 0.7)
dev.off()