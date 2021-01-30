#Use Poisson regression remove systematic biases in Hi-C trans contact maps
#Ming Hu (minghu@fas.harvard.edu)
#Last update: 08.05.2012

#read in input file
par(cex=4)
require('pheatmap')
for(chr1 in 1:23){
  for(chr2 in (chr1+1):24){

u1<-read.table(paste0("./matrix/trans/raw/chr",
                      chr1,"_",chr2,"_200kb_rawmatrix.txt"))      #user can change the name of this input file
v_a1<-read.table(paste0("/lustre/user/liclab/lirf/Project/hic/genomeFeature/mobl/200000/chr",
                       chr1,".txt",collapse=""),head=T)  #user can change the name of this input file
v_b1<-read.table(paste0("/lustre/user/liclab/lirf/Project/hic/genomeFeature/mobl/200000/chr",
                       chr2,".txt",collapse=""),head=T)  #user can change the name of this input file
a1 = intersect(intersect(which(v_a1[,4]!=0),which(v_a1[,5]!=0)),which(v_a1[,6]!=0))
b1 = intersect(intersect(which(v_b1[,4]!=0),which(v_b1[,5]!=0)),which(v_b1[,6]!=0))
v_a = v_a1[a1,]
v_b = v_b1[b1,]
u=u1[a1,b1]
#change matrix into vector
u_vec<-c(as.matrix(u))

#get cov matrix
len_m<-as.matrix(log(v_a[,4]%o%v_b[,4]))
gcc_m<-as.matrix(log(v_a[,5]%o%v_b[,5]))
map_m<-as.matrix(log(v_a[,6]%o%v_b[,6]))

#centralize cov matrix of enz, gcc
len_m<-(len_m-mean(c(len_m)))/sd(c(len_m))
gcc_m<-(gcc_m-mean(c(gcc_m)))/sd(c(gcc_m))

#change matrix into vector
len_vec<-c(len_m)
gcc_vec<-c(gcc_m)
map_vec<-c(map_m)

#fit Poisson regression: u~len+gcc+offset(map)
fit<-glm(u_vec~len_vec+gcc_vec+offset(map_vec),family="poisson")

#user can use the following two lines to fit negative binomial regression: u~len+gcc+offset(map).
#library("MASS")
#fit<-glm.nb(u_vec~len_vec+gcc_vec+offset(map_vec))

#summary(fit)
coeff<-round(fit$coeff,4)
res<-round(u/exp(coeff[1]+coeff[2]*len_m+coeff[3]*gcc_m+map_m), 4)
nor_matrix = matrix(nrow=dim(u1)[1],ncol=dim(u1)[2])
nor_matrix[a1,b1] <- as.matrix(res)
colnames(nor_matrix) = v_b1[,2]
rownames(nor_matrix) = v_a1[,2]

#output normalized trans contact map, user can change the name of this output file
write.table(res, file=paste0("./matrix/trans/normalization/chr",chr1,"_",chr2,"_normalized_matrix.txt",collapse=""), 
            row.names=F, col.names=F, sep="\t", quote=F)
pdf(paste0("./Graph/heatmapForRawMarix/trans/chr",chr1,"_",chr2,"_raw.pdf",collapse=""),onefile = F)
pheatmap(as.matrix(log(u1+1)),main = paste0("Heatmap.raw.chr",chr1,"_",chr2,collapse=""),fontsize_row = 2,cluster_rows = F, cluster_cols = F,
         border_color = NA,fontsize_col = 2,color = colorRampPalette(c( "white","red"))(50))
dev.off()
pdf(paste0("./Graph/heatmapForNormalizedMarix/trans/chr",chr1,"_",chr2,"_normalized.pdf",collapse=""),onefile = F)
pheatmap(as.matrix(log(nor_matrix+1)),main = paste0("Heatmap.normalized.chr",chr1,"_",chr2,collapse=""),fontsize_row = 2,cluster_rows = F, cluster_cols = F,
         border_color = NA,fontsize_col = 2,color = colorRampPalette(c( "white","red"))(50))
dev.off()

  }
}


