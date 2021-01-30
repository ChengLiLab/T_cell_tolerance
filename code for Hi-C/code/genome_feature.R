# 计算GC content
setwd("/lustre/user/liclab/lirf/Project/hic/data.2015.6.24/release10/normalization/mobi/tmp.data/mappabilityForhindIII")
options(scipen=3)
library(BSgenome)
library(Biostrings)
installed.genomes()
library(BSgenome.Hsapiens.UCSC.hg19)
hg19 <- BSgenome.Hsapiens.UCSC.hg19

library(GenomicRanges)

hind3 <- 'GATC'

resolution = 200000

for(chrnum in 1:24){
bin = as.integer(length(hg19[[chrnum]])/resolution,0)+1
# 找出hind3 在chr1上所有的酶切位点
hind3.chr1 <- matchPattern(hind3, hg19[[chrnum]])
hind3.start <- start(hind3.chr1)


gr <- GRanges(seqnames = 'chr1', strand = '+',
              ranges = IRanges(start = (start(hind3.chr1) - bin), width = 500))
gr.reduce <- reduce(gr)

gr.reduce[start(gr.reduce) < resolution]


## gc contant:
gc.ratio <- rep(0, bin)
for (i in c(0:(bin-1))){
  gr.seq <- gr.reduce[i * resolution  < start(gr.reduce) & start(gr.reduce) < (i + 1) * resolution ]
  seq <-  Views(hg19[[chrnum]], start = start(gr.seq), width = 400) 
  seq.gc <- alphabetFrequency(seq, baseOnly =T)
  gc.ratio[i] <- (sum(seq.gc[, 2]) + sum(seq.gc[, 3])) / sum(seq.gc)
}

## effect length:
effect.len <- rep(0, bin)
gr2 <- GRanges(seqnames = 'chr1', strand = '+',
              ranges = IRanges(start = (start(hind3.chr1) - 500), width = 1000))
gr2.reduce <- reduce(gr2)
width(gr2.reduce)

effect.len <- rep(0, bin)

for (i in c(0:(bin-1))){
  
  gr2.seq <- gr2.reduce[i * resolution  < start(gr2.reduce) & start(gr2.reduce) < (i + 1) * resolution ]
 
  effect.len[i] <- sum( width(gr2.seq) ) 
}

map = read.table(paste0("mappability.chr",chrnum,".txt",collapse=""))
feature = cbind(rep(chrnum,bin),seq(0,by=resolution,length=bin),seq(resolution,by=resolution,length=bin),effect.len,gc.ratio,map)
colnames(feature) = c("chr","start","end","len","gcc","map")
write.table(feature,file=paste0("/lustre/user/liclab/lirf/Project/hic/genomeFeature/mobl/",
                                resolution,"/chr",chrnum,".txt",collapse=""),row.names=F)
}



