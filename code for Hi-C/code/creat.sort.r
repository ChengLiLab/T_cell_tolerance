#!/usr/bin/env Rscript
# 此脚本用于通过SAMTOOLs工具将BAM文件进行排序得到sorted BAM file

files_dir <- "./tmp.data/BAM" # 要排序的BAM文件位置

files_path <- dir(files_dir, pattern = ".bam", full.names = TRUE)

n <- length(files_path)

sink("./tmp.code/sortbam.sh")

cat("#!/usr/bin/env bash\n\n")

for (i in 1:n){
cmd <- sprintf("samtools sort  -m 2000000000  %s -o ./tmp.data/BAM/%s  &",files_path[i],paste0(gsub(".bam", "", basename(files_path[i])), "sorted.bam.bam"))
# samtools sort : 将bam文件排序
# -m INT: 大约用到的最大内存，此处为2G
cat(cmd, "\n")

if (i %% 10 == 0){
            cat("wait\n")
	        }
}
cat("wait\n")
cat("sleep 3\n")

sink()
