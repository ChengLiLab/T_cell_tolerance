#!/usr/bin/env Rscript
# 将经过unique, sort ，rmp ，merge之后的bam文件转化为sam文件
files_dir <- "./tmp.data/BAM"
files_path <- dir(files_dir, pattern = "namesorted.bam.bam$", full.names = TRUE)
n <- length(files_path)
sink("./tmp.code/bamtosam.sh")
cat("#!/usr/bin/env bash\n\n")
for (i in 1:n){
cmd <- sprintf("samtools view    %s > %s  &",files_path[i],paste0("./tmp.data/BAM/",gsub("samrmp.namesorted.bam.bam", "", basename(files_path[i])), "sam"))
cat(cmd, "\n")
    if (i %% 5 == 0){
            cat("wait\n")
	        }
}
cat("wait\n")
cat("sleep 3\n")
sink()
