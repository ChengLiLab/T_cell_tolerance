#!/usr/bin/env Rscript

GTF <- "/lustre/user/liclab/publicData/igenomes/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf"
INDEX <- "/lustre/user/liclab/publicData/igenomes/Homo_sapiens/Ensembl/GRCh37/Sequence/Bowtie2Index/genome"

op.dir <- file.path(getwd(), "Tophat.out")
fq.path <- file.path(getwd())
fq.path <- file.path(fq.path, dir(fq.path, pattern = ".fastq.gz|.fq.gz"))

options(stringsAsFactors = FALSE)

sink("cufflinks.run.sh")
cat("#!/usr/bin/env bash\n\n\n")

ind <- unique(gsub("*gz$", "", fq.path))

for (i in 1:length(ind)){
if(i %% 2 != 0){
	op <- file.path(op.dir, basename(ind[i]))
	r <- fq.path[grepl(ind[i], fq.path)]
	m <- fq.path[grepl(ind[i+1], fq.path)]
	cmd <- sprintf("cufflinks -p 5 -G %s -o %s/cufflinks %s/accepted_hits.bam   &", GTF, op, op)
	cat(cmd, "\n")
	if ((i+1) %% 12 == 0){
	cat("wait\n")
	}
}}
sink()
