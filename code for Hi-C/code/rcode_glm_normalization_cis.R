#Use Poisson regression remove systematic biases in Hi-C cis contact maps
#Ming Hu (minghu@fas.harvard.edu)
#Last update: 08.05.2012
#Last update: 07.07.2015 by Ruifeng Li
#read in input file
par(cex=4)
require('pheatmap')

for(chr in 1:24){

u1<-read.table(paste0("./matrix/cis/raw/chr",
                      chr,"_200kb_rawmatrix.txt",collapse=""))           #user can change the name of this input file
v1<-read.table(paste0("/lustre/user/liclab/lirf/Project/hic/genomeFeature/mobl/200000/chr",
                      chr,".txt",collapse=""),head=T) #user can change the name of this input file
u = u1[which(v1[,4]!=0),which(v1[,4]!=0)]
v = v1[which(v1[,4]!=0),]
u = u[which(v[,5]!=0),which(v[,5]!=0)]
v = v[which(v[,5]!=0),]
u = u[which(v[,6]!=0),which(v[,6]!=0)]
v = v[which(v[,6]!=0),]
v = as.matrix(v)
bin = dim(u1)[1]
for(i in 1:dim(v)[1]){u[i,i]=0}

#change matrix into vector
u<-as.matrix(u)
u_vec<-u[upper.tri(u,diag=F)]

#get cov matrix
len_m<-as.matrix(log(v[,4]%o%v[,4]))
gcc_m<-as.matrix(log(v[,5]%o%v[,5]))
map_m<-as.matrix(log(v[,6]%o%v[,6]))

#centralize cov matrix of enz, gcc
len_m<-(len_m-mean(c(len_m)))/sd(c(len_m))
gcc_m<-(gcc_m-mean(c(gcc_m)))/sd(c(gcc_m))

#change matrix into vector
len_vec<-len_m[upper.tri(len_m,diag=F)]
gcc_vec<-gcc_m[upper.tri(gcc_m,diag=F)]
map_vec<-map_m[upper.tri(map_m,diag=F)]

#fit Poisson regression: u~len+gcc+offset(map)
fit<-glm(u_vec~len_vec+gcc_vec+offset(map_vec),family="poisson")

#user can use the following two lines to fit negative binomial regression: u~len+gcc+offset(map).
#library("MASS")
#fit<-glm.nb(u_vec~len_vec+gcc_vec+offset(map_vec))

summary(fit)
coeff<-round(fit$coeff,4)
res<- round(u/exp(coeff[1]+coeff[2]*len_m+coeff[3]*gcc_m+map_m), 4)
nor_matrix = array(dim=dim(u1))
cm = intersect(intersect(which(v1[,4]!=0),which(v1[,5]!=0)),which(v1[,6]!=0))
nor_matrix[cm,cm] = res
colnames(nor_matrix) = v1[,2]
rownames(nor_matrix) = v1[,2]
#output normalized cis contact map, user can change the name of this output file
write.table(nor_matrix, file=paste0("./matrix/cis/normalization/chr",chr,"_normalized_matrix.txt",collapse=""), 
            row.names=F, col.names=F, sep="\t", quote=F)
pdf(paste0("./Graph/heatmapForRawMarix/chr",chr,"_raw.pdf",collapse=""),width = 10,height = 10,onefile = F)
pheatmap(as.matrix(log(u1+1)),main = paste0("Heatmap.raw.chr",chr,collapse=""),fontsize_row = 2,cluster_rows = F, cluster_cols = F,
         border_color = NA,fontsize_col = 2,color = colorRampPalette(c( "white","red"))(50))
dev.off()
pdf(paste0("./Graph/heatmapForNormalizedMarix/chr",chr,"_normalized.pdf",collapse=""),width = 10,height = 10,onefile = F)
pheatmap(as.matrix(log(nor_matrix+1)),main = paste0("Heatmap.normalized.chr",chr,collapse=""),fontsize_row = 2,cluster_rows = F, cluster_cols = F,
         border_color = NA,fontsize_col = 2,color = colorRampPalette(c( "white","red"))(50))
dev.off()

}