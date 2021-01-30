library(diffHic)
library(reshape2)
library(Homo.sapiens)
human.genes = genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
seg.frags = tileGenome(seqinfo(human.genes), tilewidth=40000, cut.last.tile.in.chrom=T)
system("mkdir TAD_files")
resolution_ = "40kb"

##########################
for(num in c(1:23)){
  if(num == 23){num = "X"}
  finder =  seg.frags[seqnames(seg.frags) == paste0("chr",num,collapse = "")]
  
  name = lapply("num", paste,"|hg19|",as.character(seqnames(finder))[1],":",start(finder),"-",end(finder)+1,sep="")[[1]]
  name1 = name
  name1[1] = paste0("\t",name1[1],collapse = "")
  if(num == "X"){num = 23}
  Matrix = read.table(paste(c("./ice_normalization/chr",num,"_40k_normalized_matrix.txt"),collapse=''))

  rownames(Matrix) = name
  
  colnames(Matrix) = name1
  
  write.table(Matrix,file = paste0("./TAD_files/chr",num,"_chr",num,"_",resolution_,"_normalmatrix.txt",collapse = ""),quote = F,sep = "\t")
  
  print(num)
}





