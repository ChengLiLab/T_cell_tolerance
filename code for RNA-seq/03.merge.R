path =  dir()

gene_bio = read.table("/lustre/user/liclab/lirf/Project/sunyj/RNA.and.CHIP.seq.data/2015.6.15/data/Tophat.out/Sample_PR10_SYJ_GFP-1.1.clean.fq.111/cufflinks/genes.fpkm_tracking",header = T) 
gene_exp = read.table(paste0("./",path[1],"/cufflinks/genes.fpkm_tracking",collapse = ""),header = T)
gene = intersect(as.character(unique(gene_bio[,4])),as.character(gene_exp[,5]))
gene_exp = gene_exp[match(gene,as.character(gene_exp[,5])),c(7,10)]
for(i in 2:length(path)){
  gene_exp1 = read.table(paste0("./",path[i],"/cufflinks/genes.fpkm_tracking",collapse = ""),header = T,fill = T)
  gene_exp = cbind(gene_exp,gene_exp1[match(gene,as.character(gene_exp1[,5])),10])
  print(i)
}

colnames(gene_exp) = c("loci", gsub("_R1.fq.","",gsub("Sample_PR10_SYJ_","",path[1:length(path)])))
rownames(gene_exp) = gene
#colnames(gene_exp ) = c("loci","65N","65T")


write.table(gene_exp,file = "gene_expression1.txt")
