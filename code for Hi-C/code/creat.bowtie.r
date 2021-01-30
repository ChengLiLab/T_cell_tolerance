#!/usr/bin/env Rscript
# tested by litt 20150528
# bowtie  <- "/lustre/user/dengl/tools/bowtie2-2.2.4/bowtie2"

# 设置参数：
fqfiles <-  read.table("./parameters.txt")

# fqfiles # 1: 比对基因组index, 注意最后要写到文件名前缀； 
          #2： 要比对的fastq数据的存储路径； 3： 。。备用
        # 4： 染色体1； 5： 染色体2

files_dir <- as.character(fqfiles[2,]) # 要比对的fastq数据的存储路径

#files_path <- dir(files_dir, pattern = ".trim.fq$", full.names = TRUE) 
files_path <- dir(".", pattern = ".trim.fq$", full.names = TRUE) 
#dir() 函数： 用来生成指定目录下，文件的路径， path: 路径名，默认当前路径；pattern: 正则表达式，只抓取满足正则表达式的文件名
# full.names：如果真，返回完整路径，如果假，只返回文件名

n <- length(files_path) # 要比对的文件个数

sink("./tmp.code/bowtie2.sh") # 将输出导入到某个链接或者文件，常用cat, print函数

cat("#!/usr/bin/env bash\n\n")
for (i in 1:n){
cmd <- sprintf("bowtie2 -p 6 --phred33 -x %s  %s -S ./tmp.data/%s > ./tmp.data/%s 2>&1 &",
               fqfiles[1,],files_path[i],paste0(basename(files_path[i]),".sam"), paste0( basename(files_path[i]),".log" ))
# 此处用来生成用bowtie2做序列比对的批处理sh文件
# 参数意义： -p: parallel, 平行计算，用几个线程来计算
#            --phred33: 测序数据质量标准，Input qualities are ASCII chars equal to the Phred quality plus 33.
#            -x： index, The basename of the index for the reference genome. 是比对索引的basename
#            -S：samfile, 要输出产生的sam文件

cat(cmd, "\n")
    if (i %% 3 == 0){
            cat("wait\n") # wait 命令：前面命令结束之后，wait，保证每次最多比对3组
	        }

}
cat("wait\n")
cat("sleep 3\n") # sleep很重要，他保证整个脚本执行结束之后，再执行下一个脚本，如果没有这一行，则此脚本未结束就执行下一脚本。
sink()
