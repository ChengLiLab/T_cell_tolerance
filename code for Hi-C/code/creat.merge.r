#!/usr/bin/env Rscript

FQ.dir <- './tmp.data/BAM'
FQfiles <- dir(FQ.dir, pattern = "rmp.bam$", full.names = TRUE)
n <- length(FQfiles)
sink("./tmp.code/merge.sh")

        for (i in 1:n){
                    if(i %% 2 != 0){
                                        cmd <- sprintf("samtools merge ./tmp.data/BAM/%s %s %s  2>&1 &",paste0(gsub(".bam", "", basename(FQfiles[i])), ".merge.bam"),FQfiles[i],FQfiles[i+1])
                    # 
                                        cat(cmd, "\n")
                    if (i %% 9 == 0){
                                cat("wait\n")
                    }
}                                                }
cat("wait\n")
cat("sleep 3\n")
sink()


