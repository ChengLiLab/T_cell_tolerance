for i in `ls *rmdup.bam`
do
{
	samtools index $i
}&
done
