cutadapt -m 36 -q 20,15 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAA -o ./cleandata/CTCF-input_R1clean.fq.gz -p ./cleandata/CTCF-input_R2clean.fq.gz ./CTCF-input_R1.fq.gz ./CTCF-input_R2.fq.gz & 
cutadapt -m 36 -q 20,15 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAA -o ./cleandata/CTCF-Nu_R1clean.fq.gz -p ./cleandata/CTCF-Nu_R2clean.fq.gz ./CTCF-Nu_R1.fq.gz ./CTCF-Nu_R2.fq.gz & 
cutadapt -m 36 -q 20,15 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAA -o ./cleandata/STAT-Nu_R1clean.fq.gz -p ./cleandata/STAT-Nu_R2clean.fq.gz ./STAT-Nu_R1.fq.gz ./STAT-Nu_R2.fq.gz & 
cutadapt -m 36 -q 20,15 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAA -o ./cleandata/STAT3-input_R1clean.fq.gz -p ./cleandata/STAT3-input_R2clean.fq.gz ./STAT3-input_R1.fq.gz ./STAT3-input_R2.fq.gz & 
wait
sleep 3
