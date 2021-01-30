#!/usr/bin/env Rscript
# 此文本用来生成samtools ,将sam文件转化成bam文件的批处理sh文件。

# samtools : 这是一个非常强大的工具包，sam文件的处理基本上离不开
# bcftools： 与samtools在一起的一个包，用来call variants和操作VCF, BCF文件
# 这两个工具可以使我们从SAM格式的比对文件得到VCF格式的variant call 文件。

files_dir <- "./tmp.data" # 上一步用bowtie比对好的sam文件的位置

files_path <- dir(files_dir, pattern = "[fastq,fq].sam", full.names = TRUE) # 抓取要作转化的sam文件

n <- length(files_path) # 文件的数量
sink("./tmp.code/samTobam.sh") # 将下面的输出导入到此链接中去

cat("#!/usr/bin/env bash\n\n")
for (i in 1:n){
cmd <- sprintf("samtools view -bS   %s > %s  &",files_path[i],paste0("./tmp.data/BAM/",gsub(".fastq.sam", "", basename(files_path[i])), ".bam"))
# samtools view -bS: 将SAM文件转化成二进制格式的BAM文件
# view: 提取&打印比对好的sam或者bam文件
# -b: 输出文件是BAM格式
# -S： 输入文件是SAM格式， 如果@SQ头文件缺失，则需要-t选项


cat(cmd, "\n")
    if (i %% 5 == 0){
            cat("wait\n")
	        }
}
cat("wait\n")
cat("sleep 3\n")
sink()
