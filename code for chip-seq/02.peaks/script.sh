
### merge peaks
#cat ../peaks/P14_peaks.bed  ../peaks/P26_peaks.bed  ../peaks/P4_peaks.bed /lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/20170119/peaks/P5_peaks.bed /lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/20170119/peaks/P26.1_peaks.bed /lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/20170119/peaks/P36_peaks.bed | sort -k1,1 -k2,2n - | mergeBed -i - > merged_peaks.bed

### get reads per peak
multiBamCov -bams ../align/P14.sort.rmdup.bam  ../align/P26.sort.rmdup.bam  ../align/P4.sort.rmdup.bam /lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/20170119/align/P5.sort.rmdup.bam /lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/20170119/align/P26.1.sort.rmdup.bam /lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/20170119/align/P36.sort.rmdup.bam -bed merged_peaks.bed > merged_peaks.reads.bed

awk '{print $1":"$2"_"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' merged_peaks.reads.bed > merged_peaks.reads_reform.bed



multiBamCov -bams ../align/P14.sort.rmdup.bam  ../align/P26.sort.rmdup.bam  ../align/P4.sort.rmdup.bam /lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/20170119/align/P5.sort.rmdup.bam /lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/20170119/align/P26.1.sort.rmdup.bam /lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/20170119/align/P36.sort.rmdup.bam -bed HUVEC.SE.bed.txt > SE.reads.bed
awk '{print $1":"$2"_"$3"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11}' SE.reads.bed > SE.reads_reform.bed


