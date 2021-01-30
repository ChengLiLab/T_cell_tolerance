#/usr/bin/env Rscript
# 此脚本用于去除PCR产生的duplicates

files_dir <- "./tmp.data/BAM"
files_path <- dir(files_dir, pattern = "unique.bam", full.names = TRUE)
n <- length(files_path)
sink("./tmp.code/rmpbam.sh")
cat("#!/usr/bin/env bash\n\n")
for (i in 1:n){
# cmd <- sprintf("samtools rmdup -s %s  ./tmp.data/BAM/%s  &",files_path[i],paste0(gsub("sortedunique.bam", "", basename(files_path[i])), "rmp.bam"))
# rmdup: 如果多个read pairs 含有相同的比对信息，只保留质量最高的那一个
# 默认是paired reads ， 如果是单端比对，用-s参数
# samtools用来做去除duplicate没有picard广泛，故改为用picard中的MarkDpulicates函数：
cmd <- sprintf("java -Xms2g -Xmx8g -XX:ParallelGCThreads=3 -jar /lustre/user/liclab/biotools/picard-tools-1.118/MarkDuplicates.jar INPUT=%s OUTPUT=./tmp.data/BAM/%s METRICS_FILE=./tmp.data/BAM/%s_metrics VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true REMOVE_DUPLICATES=true &", 
               files_path[i], paste0(gsub("sortedunique.bam", "", basename(files_path[i])), "rmp.bam"), paste0(gsub("sortedunique.bam", "", basename(files_path[i])), "rmp.bam"))
cat(cmd, "\n")
    if (i %% 10 == 0){
            cat("wait\n")
	        }
}
cat("wait\n")
cat("sleep 3\n")
sink()
