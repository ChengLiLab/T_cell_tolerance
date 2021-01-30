#!/usr/bin/env Rscript
# 此脚本用于将bam文件根据reads名字来排序
FQ.dir <- "./tmp.data/BAM"
FQfiles <- dir(FQ.dir, pattern = "merge.bam", full.names = TRUE)
n <- length(FQfiles)
sink("./tmp.code/sort.sh")

for (i in 1:n){
cmd <- sprintf("samtools sort -n  %s -o ./tmp.data/BAM/%s  2>&1 &",FQfiles[i],paste0(gsub(".merge.bam", "", basename(FQfiles[i])), ".namesorted.bam.bam"))
                    cat(cmd, "\n")
                    if (i %% 9 == 0){
                                cat("wait\n")
                    }
}
cat("wait\n")
cat("sleep 3\n")
sink()


