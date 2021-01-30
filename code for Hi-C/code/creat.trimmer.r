#!/usr/bin/env Rscript

# 设置参数：
fqfiles <-  read.table("./parameters.txt")

files_dir <- as.character(fqfiles[2,]) # 要比对的fastq数据的存储路径

files_path <- dir(files_dir, pattern = ".gz$", full.names = TRUE) 


n <- length(files_path) # 要比对的文件个数

sink("./tmp.code/trimmer.sh") # 将输出导入到某个链接或者文件，常用cat, print函数

cat("#!/usr/bin/env bash\n\n")
for (i in 1:n){
cmd <- sprintf("zcat %s | fastx_trimmer  -Q33  -l 36  -o %s  &",
               files_path[i],paste0(gsub("clean.fq", "", basename(files_path[i])),"clean.trim.fq"))



cat(cmd, "\n")
    if (i %% 4 == 0){
            cat("wait\n") # wait 命令：前面命令结束之后，wait，保证每次最多比对3组
	        }

}
cat("wait\n")
cat("sleep 3\n") # sleep很重要，他保证整个脚本执行结束之后，再执行下一个脚本，如果没有这一行，则此脚本未结束就执行下一脚本。
sink()
