#/usr/bin/env Rscript
files_dir <- "./tmp.data/BAM"
files_path <- dir(files_dir, pattern = "^sam_merge.txt$", full.names = TRUE)
n <- length(files_path)
sink("./tmp.code/filter.py.sh")
cat("#!/usr/bin/env bash\n\n")
for (i in 1:n){
cmd <- sprintf("python /lustre/user/liclab/lirf/Project/hic/bgi/5534N/code/filter.py  -i  %s -o ./tmp.data/BAM/%s  &",files_path[i],paste0(gsub("sam_merge.txt", "", basename(files_path[i])), "filter.txt"))
cat(cmd, "\n")
    if (i %% 10 == 0){
            cat("wait\n")
	        }

}
cat("wait\n")
cat("sleep 3\n")
sink()
