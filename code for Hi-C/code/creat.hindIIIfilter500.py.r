#/usr/bin/env Rscript
files_dir <- "./tmp.data/BAM"
files_path <- dir(files_dir, pattern = "^filter.txt$", full.names = TRUE)
n <- length(files_path)
sink("./tmp.code/hindIIIfilter500.py.sh")
cat("#!/usr/bin/env bash\n\n")
for (i in 1:n){
cmd <- sprintf("python /lustre/user/liclab/lirf/Project/hic/bgi/5534N/code/hindIIIfilter500.py  -i  %s -o ./tmp.data/BAM/%s  &",files_path[i],paste0(gsub(".txt", "", basename(files_path[i])), ".hindIII.txt"))
cat(cmd, "\n")
    if (i %% 10 == 0){
            cat("wait\n")
	        }

}
cat("wait\n")
cat("sleep 3\n")
sink()
