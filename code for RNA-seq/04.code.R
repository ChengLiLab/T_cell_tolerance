setwd("/lustre/user/liclab/lirf/Project/renmin-chang/p2/rnaseq")
exp1 = read.table("/lustre/user/liclab/lirf/Project/renmin-chang/p2/before/rna/gene_expression1.txt",header = T,row.names = 1,sep = " ")
exp2 = read.table("/lustre/user/liclab/lirf/Project/renmin-chang/p2/after/rna/gene_expression1.txt",header = T,row.names = 1,sep = " ")
exp = cbind(exp1,exp2[,2:4])
# row.names(exp)[which((is.na(exp[,3])==T))]
# exp[which((is.na(exp[,8])==T)),8]
# mis = match(c("CSH2","CFL1","KIAA1432","ERMP1","ATP11A","CD274","RAPH1","C20orf194","SLC19A3","ABI2"),row.names(exp))
# mis = c(which((is.na(exp[,8])==T)),mis)
# #exp = exp[-mis,]
# exp[,3] = as.numeric(as.character(exp[,3]))
# exp[,8] = as.numeric(as.charactefr(exp[,8]))
# exp[which((is.na(exp[,8])==T)),8] = 0
# exp[,13] = as.numeric(as.character(exp[,13]))

##############boxplot

colnames(exp) = c("loci","P4-CD4_before","P4-CD8_before","P4-NK_before","P4-CD4_after","P4-CD8_after","P4-NK_after")

#exp_p1_p2 = exp
exp_p1_p2 = as.matrix(exp_p1_p2)
exp = cbind(exp_p1_p2,exp_p3[,-1],exp_p4[,-1])
colnames(exp)[22] = "P4-NK_after"
colnames(exp)
par(mar=c(10,6,4,2))
boxplot(exp2[,1:22],ylim = c(0,30),main="",ylab = "FPKM",col = c("steelblue","steelblue","steelblue",
                                                                "hotpink","hotpink","hotpink"),las = 2)
#save(exp,file = "all_exp.dat")


##########pca 
##exp = read.table("sun_expression.2016.1.11.csv",header = T,row.names = 1,sep = ";")
exp = as.matrix(exp)
sd_exp = apply(exp,1,sd)
order(sd_exp)[1:100]
###TOP10 gene
top10 = sd_exp[order(sd_exp,decreasing = T)[1:1000]]
pheatmap(as.matrix(log(exp[names(top10),]+1)),main = "heatmap",cluster_rows = T, cluster_cols = T,fontsize_row = 1,
         border_color = NA,color = colorRampPalette(c( "blue","white","red"))(50))

source("/lustre/user/liclab/lirf/Project/shuoshuo/rnaseq/normalize.r")
exp3 = quantileNormalizeByFeature(matrix_to_normalize = t(exp[,13:22]),target_distribution_matrix = t(exp[,1:12]) )
exp2 = as.matrix(cbind(t(exp3),exp[,1:12]))
sd_exp = apply(exp2,1,sd)
order(sd_exp)[1:100]
top10 = sd_exp[order(sd_exp,decreasing = T)[1:2000]]
pheatmap(as.matrix(log(exp2[names(top10),]+1)),main = "heatmap",cluster_rows = T, cluster_cols = T,fontsize_row = 1,
         border_color = NA,color = colorRampPalette(c( "blue","white","red","red"))(50))


##########pca 
##exp = read.table("sun_expression.2016.1.11.csv",header = T,row.names = 1,sep = ";")
#exp = as.matrix(exp[,2:13])
sd_exp = apply(cleandat,1,sd)
order(sd_exp)[1:100]
exp_2000 = cleandat[order(sd_exp,decreasing = T)[1:2000],]

b=matrix(as.numeric(exp_2000),nrow=nrow(exp_2000))
colnames(b) = colnames(exp_2000)
rownames(b) = rownames(exp_2000)
#####cluster
hc <- hclust(dist(t(log(b+1))))
plot(hc, hang = -1)
plclust( hc, hang = -1)
re <- rect.hclust(hc, k = 4)
iris.id <- cutree(hc, 3)

save













exp_2000 = exp[order(sd_exp,decreasing = T)[10:110],]

b=matrix(as.numeric(exp_2000),nrow=nrow(exp_2000))
colnames(b) = colnames(exp_2000)
rownames(b) = rownames(exp_2000)
#####cluster
hc <- hclust(dist(t(log(b+1))))
plot(hc, hang = -1)
plclust( hc, hang = -1)
re <- rect.hclust(hc, k = 3)
iris.id <- cutree(hc, 3)
#########基因的热图
library(pheatmap)
pheatmap(as.matrix(log(b+1)),main = "heatmap",cluster_rows = T, cluster_cols = T,fontsize_row = 6,
         border_color = NA,color = colorRampPalette(c( "black","white","orange"))(50))

























#exp = exp[,2:7]


###########overlap
cd4_rna = read.csv("/lustre/user/liclab/lirf/Project/renmin-chang/analysis_20181006/after/rnaseq/Tophat.out/cd4.csv")






result = array(NA,dim = c(dim(exp)[1],7))
colnames(exp)
result[,1] = exp[,"P4-CD4_after"]/exp[,"P4-CD4_before"]
result[,2] = abs(exp[,"P4-CD4_after"]-exp[,"P4-CD4_before"])
result[,3] = exp[,"P4-CD8_after"]/exp[,"P4-CD8_before"]
result[,4] = abs(exp[,"P4-CD8_after"]-exp[,"P4-CD8_before"])
result[,5] = exp[,"P4-NK_after"]/exp[,"P4-NK_before"]
result[,6] = abs(exp[,"P4-NK_after"]-exp[,"P4-NK_before"])
result[,7] = apply(exp[,17:22], 1, var)

rr = cbind(exp[,17:22],result)
write.table(rr,file = "result1.txt",col.names = T,row.names = T,quote = F,sep = "\t")
quantile(result[,7])

########CD4
colnames(exp)
LS = c("ZZSG.BMCD3CD4","LZDG.BMCD3CD4","P3-CD4_after","P4-CD4_after")
NL = c("LZDNG.BMCD3CD4","ZZSNG.BMCD3CD4","P3-CD4_before","P4-CD4_before")
exp[which(is.na(exp[,"LZDG.BMCD3CD4"]) == T),"LZDG.BMCD3CD4"] = 0


LS= as.character(LS)
NL= as.character(NL)

result = array(NA,dim = c(dim(exp)[1],3))
for( i in 1:dim(exp)[1]){
  result[i,1]  = mean(as.numeric(exp[i,LS])) - mean(as.numeric(exp[i,NL]))
  if(mean(as.numeric(exp[i,LS])) > 1)
  {result[i,2]  = t.test(as.numeric(exp[i,LS]), as.numeric(exp[i,NL]))$p.value}
  print(i)
}




result[,3] = p.adjust(result[,2])
rownames(result) = rownames(exp)
colnames(result) = c("CD4_FC", "CD4_P-value" , "CD4_adj.P-value")
write.csv(result,file = "CD4_result.csv")





##########pca 
##exp = read.table("sun_expression.2016.1.11.csv",header = T,row.names = 1,sep = ";")
exp = as.matrix(exp[,2:13])
sd_exp = apply(exp,1,sd)
order(sd_exp)[1:100]
###TOP10 gene
top10 = sd_exp[order(sd_exp,decreasing = T)[1:10]]
pheatmap(as.matrix(log(exp[names(top10),]+1)),main = "heatmap",cluster_rows = T, cluster_cols = F,fontsize_row = 20,
         border_color = NA,color = colorRampPalette(c( "black","white","orange"))(50))
exp_2000 = exp[order(sd_exp,decreasing = T)[10:110],]

b=matrix(as.numeric(exp_2000),nrow=nrow(exp_2000))
colnames(b) = colnames(exp_2000)
rownames(b) = rownames(exp_2000)
#####cluster
hc <- hclust(dist(t(log(b+1))))
plot(hc, hang = -1)
plclust( hc, hang = -1)
re <- rect.hclust(hc, k = 3)
iris.id <- cutree(hc, 3)
#########基因的热图
library(pheatmap)
pheatmap(as.matrix(log(b+1)),main = "heatmap",cluster_rows = T, cluster_cols = T,fontsize_row = 6,
         border_color = NA,color = colorRampPalette(c( "black","white","orange"))(50))



#############差异基因


#######CD3CD4
fc = (exp[,"LZDG.BMCD3CD4"]+0.0001)/(exp[,"LZDNG.BMCD3CD4"]+0.0001)  #######lzd
cha = abs(exp[,"LZDG.BMCD3CD4"] - exp[,"LZDNG.BMCD3CD4"])
up_cd3cd4_lzd = intersect(which(fc>2),which(cha>2))
down_cd3cd4_lzd = intersect(which(fc<0.5),which(cha>2))

fc = (exp[,"ZZSG.BMCD3CD4"]+0.0001)/(exp[,"ZZSNG.BMCD3CD4"]+0.0001)  #######zzs
cha = abs(exp[,"ZZSG.BMCD3CD4"] - exp[,"ZZSNG.BMCD3CD4"])
up_cd3cd4_zzs = intersect(which(fc>2),which(cha>2))
down_cd3cd4_zzs = intersect(which(fc<0.5),which(cha>2))

library(Vennerable)
data<-Venn(list("LZD_UP"=up_cd3cd4_lzd,
                "ZZS_UP"=up_cd3cd4_zzs))   
plot(data,doWeight=T)   ######按面积比例画
data<-Venn(list("LZD_DOWN"=down_cd3cd4_lzd,
                "ZZS_DOWN"=down_cd3cd4_zzs))   
plot(data,doWeight=T)   ######按面积比例画

intersect(up_cd3cd4_lzd,up_cd3cd4_zzs)

write.table(rownames(exp)[intersect(up_cd3cd4_lzd,up_cd3cd4_zzs)],file = "up_cd3cd4.txt",col.names = F,row.names = F,quote = F)
write.table(rownames(exp)[intersect(down_cd3cd4_lzd,down_cd3cd4_zzs)],file = "down_cd3cd4.txt",col.names = F,row.names = F,quote = F)

#######CD3CD8
fc = (exp[,"LZDG.BMCD3CD8"]+0.0001)/(exp[,"LZDNG.BMCD3CD8"]+0.0001)  #######lzd
cha = abs(exp[,"LZDG.BMCD3CD8"] - exp[,"LZDNG.BMCD3CD8"])
up_cd3cd8_lzd = intersect(which(fc>2),which(cha>2))
down_cd3cd8_lzd = intersect(which(fc<0.5),which(cha>2))

fc = (exp[,"ZZSG.BMCD3CD8"]+0.0001)/(exp[,"ZZSNG.BMCD3CD8"]+0.0001)  #######zzs
cha = abs(exp[,"ZZSG.BMCD3CD8"] - exp[,"ZZSNG.BMCD3CD8"])
up_cd3cd8_zzs = intersect(which(fc>2),which(cha>2))
down_cd3cd8_zzs = intersect(which(fc<0.5),which(cha>2))

library(Vennerable)
data<-Venn(list("LZD_UP"=up_cd3cd8_lzd,
                "ZZS_UP"=up_cd3cd8_zzs))   
plot(data,doWeight=T)   ######按面积比例画
data<-Venn(list("LZD_DOWN"=down_cd3cd8_lzd,
                "ZZS_DOWN"=down_cd3cd8_zzs))   
plot(data,doWeight=T)   ######按面积比例画

intersect(up_cd3cd8_lzd,up_cd3cd8_zzs)

write.table(rownames(exp)[intersect(up_cd3cd8_lzd,up_cd3cd8_zzs)],file = "up_cd3cd8.txt",col.names = F,row.names = F,quote = F)
write.table(rownames(exp)[intersect(down_cd3cd8_lzd,down_cd3cd8_zzs)],file = "down_cd3cd8.txt",col.names = F,row.names = F,quote = F)


#######CD19
fc = (exp[,"LZDG.BMCD19"]+0.0001)/(exp[,"LZDNG.BMCD19"]+0.0001)  #######lzd
cha = abs(exp[,"LZDG.BMCD19"] - exp[,"LZDNG.BMCD19"])
up_cd19_lzd = intersect(which(fc>2),which(cha>2))
down_cd19_lzd = intersect(which(fc<0.5),which(cha>2))

fc = (exp[,"ZZSG.BMCD19"]+0.0001)/(exp[,"ZZSNG.BMCD19"]+0.0001)  #######zzs
cha = abs(exp[,"ZZSG.BMCD19"] - exp[,"ZZSNG.BMCD19"])
up_cd19_zzs = intersect(which(fc>2),which(cha>2))
down_cd19_zzs = intersect(which(fc<0.5),which(cha>2))

library(Vennerable)
data<-Venn(list("LZD_UP"=up_cd19_lzd,
                "ZZS_UP"=up_cd19_zzs))   
plot(data,doWeight=T)   ######按面积比例画
data<-Venn(list("LZD_DOWN"=down_cd19_lzd,
                "ZZS_DOWN"=down_cd19_zzs))   
plot(data,doWeight=T)   ######按面积比例画

intersect(up_cd19_lzd,up_cd19_zzs)

write.table(rownames(exp)[intersect(up_cd19_lzd,up_cd19_zzs)],file = "up_cd19.txt",col.names = F,row.names = F,quote = F)
write.table(rownames(exp)[intersect(down_cd19_lzd,down_cd19_zzs)],file = "down_cd19.txt",col.names = F,row.names = F,quote = F)





#########Tcell 
tc = intersect(intersect(up_cd3cd4_lzd,up_cd3cd4_zzs),intersect(up_cd3cd8_lzd,up_cd3cd8_zzs))
pheatmap(as.matrix(log(exp[tc,]+1)),main = "heatmap",cluster_rows = T, cluster_cols = T,fontsize_row = 10,
         border_color = NA,color = colorRampPalette(c( "black","white","orange"))(50))

tc = intersect(intersect(down_cd3cd4_lzd,down_cd3cd4_zzs),intersect(down_cd3cd8_lzd,down_cd3cd8_zzs))
pheatmap(as.matrix(log(exp[tc,c(1,2,3,7,8,9,4,5,6,10,11,12)]+1)),main = "heatmap",cluster_rows = T, cluster_cols = F,fontsize_row = 10,
         border_color = NA,color = colorRampPalette(c( "black","white","orange"))(50))

tc = intersect(up_cd19_lzd,up_cd19_zzs)
pheatmap(as.matrix(log(exp[tc,]+1)),main = "heatmap",cluster_rows = T, cluster_cols = T,fontsize_row = 4,
         border_color = NA,color = colorRampPalette(c( "black","white","orange"))(50))

tc = intersect(down_cd19_lzd,down_cd19_zzs)
pheatmap(as.matrix(log(exp[tc,c(1,2,3,7,8,9,4,5,6,10,11,12)]+1)),main = "heatmap",cluster_rows = T, cluster_cols = F,fontsize_row = 4,
         border_color = NA,color = colorRampPalette(c( "black","white","orange"))(50))
colnames(exp)






