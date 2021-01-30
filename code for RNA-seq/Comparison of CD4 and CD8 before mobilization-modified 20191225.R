setwd("/lustre/user/liclab/lirf/Project/renmin-chang/p3/rnaseq/动员前CD4和CD8比较")
load("/lustre/user/liclab/lirf/Project/renmin-chang/p3/rnaseq/exp.data")
########################CD4
colnames(exp)
colnames(exp) = c("HD2-CD4-Tss", "HD2-CD8-Tss", "HD2-CD4-Ttol",  "HD2-CD8-Ttol",  "HD3-CD4-Ttol",  "HD3-CD8-Ttol",  "HD3-CD4-Tss",
  "HD3-CD8-Tss", "HD1-CD4-Tss", "HD1-CD8-Tss", "HD1-CD4-Ttol","HD1-CD8-Ttol")
dim(exp)
exp[1:5,1:5]
exp_log2 = log2(exp+1)
NL = c("HD1-CD8-Tss","HD2-CD8-Tss","HD3-CD8-Tss")
LS = c("HD1-CD4-Tss","HD2-CD4-Tss","HD3-CD4-Tss")
LS= as.character(LS)
NL= as.character(NL)
LS
NL
#exp_log2
result = array(NA,dim = c(dim(exp_log2)[1],6))
for( i in 1:dim(exp_log2)[1]){
  result[i,1]  = mean(as.numeric(exp_log2[i,LS])) - mean(as.numeric(exp_log2[i,NL]))
  if(mean(as.numeric(exp_log2[i,LS])) != as.numeric(exp_log2[i,NL])){result[i,2]  = t.test(as.numeric(exp_log2[i,LS]), as.numeric(exp_log2[i,NL]))$p.value}
  result[i,3] = mean(as.numeric(exp_log2[i,NL]))
  result[i,4] = mean(as.numeric(exp_log2[i,LS]))
  result[i,5] = mean(as.numeric(exp_log2[i,c(LS,NL)]))
  print(i)
  
}
result[,6] = p.adjust(result[,2])
rownames(result) = rownames(exp_log2)
colnames(result) = c("FC(log2)", "P-value" ,"CD4(mean)","CD8(mean)", "ALL(mean)","adj.P-value")

#write.csv(result,file = "CD8vsCD4.diff.csv")
########################
colnames(result)
summary(result[,5])
result1 = result[which(result[,5] >=1),]
# top = result1[order(result1[,2],decreasing = F)[1:20],]
# top = top[order(top[,1],decreasing = F)[1:20],]
top = result1[which(result1[,2] < 0.05),]
top = top[c(order(top[,1],decreasing = F)[1:10],order(top[,1],decreasing = T)[1:10]),]
top = top[order(top[,1],decreasing = T),]

topgene = rownames(top)
aa = as.matrix(exp[topgene,c(NL,LS)])
bb = t(scale(t(aa),center = F))
#bb = scale(aa)
annotation_row = data.frame(FC_log2 = top[,1],P.value_log10 = -log10(top[,2]))
rownames(annotation_row) = rownames(bb)
ann_colors = list(FC_log2 = c("darkgreen","white", "orange"),P.value_log10 = c("white", "hotpink"))
pheatmap(bb,main = "",cluster_rows = F, cluster_cols = F,fontsize_row = 8,annotation_row = annotation_row,annotation_colors = ann_colors,
         border_color = NA,color = colorRampPalette(c( "royalblue1","white","brown1"))(50))


################################################################################动员后

NL = c("HD1-CD8-Ttol","HD2-CD8-Ttol","HD3-CD8-Ttol")
LS = c("HD1-CD4-Ttol","HD2-CD4-Ttol","HD3-CD4-Ttol")
LS= as.character(LS)
NL= as.character(NL)
LS
NL
#exp_log2
result = array(NA,dim = c(dim(exp_log2)[1],6))
for( i in 1:dim(exp_log2)[1]){
  result[i,1]  = mean(as.numeric(exp_log2[i,LS])) - mean(as.numeric(exp_log2[i,NL]))
  if(mean(as.numeric(exp_log2[i,LS])) != as.numeric(exp_log2[i,NL])){result[i,2]  = t.test(as.numeric(exp_log2[i,LS]), as.numeric(exp_log2[i,NL]))$p.value}
  result[i,3] = mean(as.numeric(exp_log2[i,NL]))
  result[i,4] = mean(as.numeric(exp_log2[i,LS]))
  result[i,5] = mean(as.numeric(exp_log2[i,c(LS,NL)]))
  print(i)
  
}
result[,6] = p.adjust(result[,2])
rownames(result) = rownames(exp_log2)
colnames(result) = c("FC(log2)", "P-value" ,"CD4(mean)","CD8(mean)", "ALL(mean)","adj.P-value")

#write.csv(result,file = "CD8vsCD4.diff.csv")
########################
colnames(result)
summary(result[,5])
result1 = result[which(result[,5] >=1),]
# top = result1[order(result1[,2],decreasing = F)[1:20],]
# top = top[order(top[,1],decreasing = F)[1:20],]
top = result1[which(result1[,2] < 0.05),]
top = top[c(order(top[,1],decreasing = F)[1:10],order(top[,1],decreasing = T)[1:10]),]
top = top[order(top[,1],decreasing = T),]

topgene = rownames(top)
aa = as.matrix(exp[topgene,c(NL,LS)])
bb = t(scale(t(aa),center = F))
#bb = scale(aa)
library(pheatmap)
annotation_row = data.frame(FC_log2 = top[,1],P.value_log10 = -log10(top[,2]))
rownames(annotation_row) = rownames(bb)
ann_colors = list(FC_log2 = c("darkgreen","white", "orange"),P.value_log10 = c("white", "hotpink"))
pheatmap(bb,main = "",cluster_rows = F, cluster_cols = F,fontsize_row = 8,annotation_row = annotation_row,annotation_colors = ann_colors,
         border_color = NA,color = colorRampPalette(c( "royalblue1","white","brown1"))(50))
########火山图



result = as.data.frame(result)
result = na.omit(result)
result$`FC(log2)`
result$`P-value`
Threshold <- as.factor((result$`FC(log2)` > 1.2 | result$`FC(log2)` < -1.2)& result$`P-value` < 0.05)
ggplot(result,aes(x = `FC(log2)`, y = -log10(`P-value`),colour = Threshold)) + 
  xlab("log2(Fold Change)")+ylab("-log10(p-value)") + geom_point()+ geom_text(label=paste(rownames(result)),colour="black",size=1)



cleandat = log2(exp+1)
sd_exp = apply(cleandat,1,sd)
exp_2000 = cleandat[order(sd_exp,decreasing = T)[1:1000],]

b=matrix(as.numeric(as.matrix(exp_2000)),nrow=nrow(exp_2000))
colnames(b) = colnames(exp_2000)
rownames(b) = rownames(exp_2000)
#####cluster
hc <- hclust(dist(t(b)))
plot(hc, hang = -1,lwd = 2)


#####################plot cluster 20201208
colnames(exp)
exp_log2 = log2(exp+1)
NL = c("HD1-CD4-Tss","HD2-CD4-Tss","HD3-CD4-Tss","HD1-CD8-Tss","HD2-CD8-Tss","HD3-CD8-Tss")
LS = c("HD1-CD4-Ttol","HD2-CD4-Ttol","HD3-CD4-Ttol","HD1-CD8-Ttol","HD2-CD8-Ttol","HD3-CD8-Ttol")
LS= as.character(LS)
NL= as.character(NL)
LS
NL
#exp_log2
result = array(NA,dim = c(dim(exp_log2)[1],6))
for( i in 1:dim(exp_log2)[1]){
  result[i,1]  = mean(as.numeric(exp_log2[i,LS])) - mean(as.numeric(exp_log2[i,NL]))
  if(mean(as.numeric(exp_log2[i,LS])) != as.numeric(exp_log2[i,NL])){result[i,2]  = t.test(as.numeric(exp_log2[i,LS]), as.numeric(exp_log2[i,NL]))$p.value}
  result[i,3] = mean(as.numeric(exp_log2[i,NL]))
  result[i,4] = mean(as.numeric(exp_log2[i,LS]))
  result[i,5] = mean(as.numeric(exp_log2[i,c(LS,NL)]))
  print(i)
  
}
result[,6] = p.adjust(result[,2])
rownames(result) = rownames(exp_log2)
colnames(result) = c("FC(log2)", "P-value" ,"CD4(mean)","CD8(mean)", "ALL(mean)","adj.P-value")




result1 = result[which(result[,5] >=1),]
# top = result1[order(result1[,2],decreasing = F)[1:20],]
# top = top[order(top[,1],decreasing = F)[1:20],]
top = result1[which(result1[,2] < 1),]
top = top[c(order(top[,1],decreasing = F)[1:440],order(top[,1],decreasing = T)[1:440]),]
top = top[order(top[,1],decreasing = T),]

topgene = rownames(top)

cleandat = log2(exp+1)
exp_2000 = cleandat[topgene,]

b=matrix(as.numeric(as.matrix(exp_2000)),nrow=nrow(exp_2000))
colnames(b) = colnames(exp_2000)
rownames(b) = rownames(exp_2000)
#####cluster
hc <- hclust(dist(t(b)))
plot(hc, hang = -1,lwd = 2)

