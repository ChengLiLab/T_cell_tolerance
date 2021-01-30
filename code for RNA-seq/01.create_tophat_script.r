#!/usr/bin/env Rscript

GTF <- "/lustre/user/liclab/publicData/igenomes/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf"
INDEX <- "/lustre/user/liclab/publicData/igenomes/Homo_sapiens/Ensembl/GRCh37/Sequence/Bowtie2Index/genome"

op.dir <- file.path(getwd(), "Tophat.out")
if (!file.exists(op.dir)) dir.create(op.dir)

fq.path <- file.path(getwd())
#fq.path <- file.path(fq.path, dir(fq.path, pattern = "gz"))
fq.path <- file.path(fq.path, dir(fq.path, pattern = ".fastq.gz|.fq.gz"))

options(stringsAsFactors = FALSE)
#si <- read.table("./sampleInfo.csv", sep = ",")

sink("tophat.run.sh")

ind <- unique(gsub("*gz$", "", fq.path))
#for (i in 1:nrow(si)){
for (i in 1:length(ind)){
    #op <- file.path(op.dir, si[i, 1])
    #r1 <- file.path(fq.path, si[i, 3])
    #r2 <- file.path(fq.path, si[i, 4])
 if(i %% 2 != 0){
    op <- file.path(op.dir, basename(ind[i]))
    r <- fq.path[grepl(ind[i], fq.path)]
    m <- fq.path[grepl(ind[i+1], fq.path)]
    cmd <- sprintf("tophat -p 4 -G %s -o %s --no-novel-juncs %s %s %s&", GTF, op, INDEX, r[1], m[1])
    cat(cmd, "\n")
    if ((i+1) %% 20 == 0){
        cat("wait\n")
    }
}}
sink()
