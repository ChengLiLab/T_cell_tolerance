inward= read.table("./tmp.data/BAM/inward.txt",header =F)
outward= read.table("./tmp.data/BAM/outward.txt",header =F)
samestrand= read.table("./tmp.data/BAM/samestrand.txt",header =F)


inw = array(dim=c(9,3))
k=1
for(i in c(0,1,10,50,100,500,1000,5000,10000)){
  inw[k,1] = length(which(inward[,1]<(100+100*i) & inward[,1]>=100*i))
  inw[k,2] = length(which(outward[,1]<(100+100*i) & outward[,1]>=100*i))
  inw[k,3] = length(which(samestrand[,1]<(100+100*i) & samestrand[,1]>=100*i))
  k=k+1
print(i)
}
par(cex=2)
pdf("./Graph/Error structure by distance.pdf",4,4)
plot(inw[,1]/inw[,3],type="l",ylim=c(0,2.5),lwd=2,axes = FALSE,
     main="Error structure by distance",xlab="Distance (gap size)",
     ylab="Ratio of reads number")
axis(2,cex = 0.7)
axis(1, seq(1,9),rep("",times = 9))
text(seq(1,9),-0.6,c("0","100bp","1kb","5kb","10kb","50kb","100kb","500kb","1Mb"),srt=-45,xpd = TRUE,cex = 0.8)
lines(inw[,2]/inw[,3],type="l",col ="red",lwd=2)
legend("topright", legend=c("inward/samestrand","outward/samestrand"), col=1:2, lwd=2,border ="white",cex = 0.7)
abline(0.5,0,lty = 2,lwd=2)
#text(7,0.25,"U226,HindIII",cex = 0.8)
dev.off()
