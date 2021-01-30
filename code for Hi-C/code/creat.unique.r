#/usr/bin/env Rscript
# 此脚本用于过滤BAM文件中的alignment
# 采用MAPQ >1 的过滤方式

files_dir <- "./tmp.data/BAM"
files_path <- dir(files_dir, pattern = ".bam.bam", full.names = TRUE)
n <- length(files_path)
sink("./tmp.code/uniquebam.sh")
cat("#!/usr/bin/env bash\n\n")
for (i in 1:n){
cmd <- sprintf("samtools view -bq 1  %s >  ./tmp.data/BAM/%s  &",files_path[i],paste0(gsub(".bam.bam", "", basename(files_path[i])), "unique.bam"))
# -b: 输出BAM格式的文件
# -q: 略过MAPQ小于0的比对， minMapQ ,此处表示 MAPQ > 1
# 不同的比对软件得到的MAPQ不同，所以最通用的过滤reads的方法是用FLAG信息来过滤，参数选择 -f或者-F
# view -f : 得到
# view -F: 去掉不满足条件的序列 ,建议BOWTIEW 在比对的时候用-k参数
cat(cmd, "\n")
    if (i %% 10 == 0){
            cat("wait\n")
	        }
}
cat("wait\n")
cat("sleep 3\n")

sink()
