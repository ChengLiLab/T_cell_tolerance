options(scipen=3)
system("mkdir ice_normalization")
source("/lustre/user/liclab/lirf/Project/hic/2014/ice_norm.R")
chr.length = read.table("/lustre/user/liclab/lirf/Project/hic/2014/25/ice/chr.length.txt")
dd = dir(path = "./raw/")
for(i in c(1:length(dd))){
u1<-read.table(paste0("./raw/",
                      dd[i],collapse=""))
chr = as.numeric(gsub("_20k_rawmatrix.txt","",gsub("chr","",dd[i])))

ice_norm = IceNorm(u1,chr,chr.length[chr,2],20000,200)
write.table(ice_norm, file=paste0("./ice_normalization/",
                                  gsub("_rawmatrix.txt","",dd[i]),"_normalized_matrix.txt",collapse=""), 
            row.names=F, col.names=F, sep="\t", quote=F)
}
