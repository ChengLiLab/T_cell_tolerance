for i in `ls *sort.rmdup.bam`
do
	bamCoverage -b $i -o ${i/sort.rmdup.bam/bw} --normalizeUsing RPKM &
done


# for i in `ls *_split/shifted.bam`
# do
# 	bamCoverage -b $i -o ${i/bam/bw} --normalizeUsingRPKM --binSize 1 --numberOfProcessors 5 &
# done

