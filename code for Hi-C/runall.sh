#!/usr/bin/env bash
###################################################################################
# pipeline for hi-c data
# 使用方法： 
# 1.将前置信息输入到parameter.txt文件中
# 2.将此文件，parameter.txt，code文件夹三个存入一个文件中，复制到要分析的数据那里
# 3. 执行此pipeline文件即可(sh runall.sh &)
###################################################################################

# 1. 序列比对
# alignment and bamfilter
mkdir ./tmp.data
mkdir ./tmp.code
mkdir ./tmp.data/BAM

# #trimmer 36bp
Rscript /lustre/user/liclab/lirf/Project/hic/bgi/5534N/code/creat.trimmer.r
sh ./tmp.code/trimmer.sh



#此R脚本用于生成bowti2比对的批处理文件 
Rscript /lustre/user/liclab/lirf/Project/hic/bgi/5534N/code/creat.bowtie.r
sh ./tmp.code/bowtie2.sh
# 
# 此脚本用于将上一步比对产生的SAM文件转化成为二进制的BAM文件
Rscript /lustre/user/liclab/lirf/Project/hic/bgi/5534N/code/creat.samtobam.r
sh ./tmp.code/samTobam.sh

# 此脚本用于将上一步产生的BAM文件排序
Rscript /lustre/user/liclab/lirf/Project/hic/bgi/5534N/code/creat.sort.r
sh ./tmp.code/sortbam.sh

# 此脚本用于产生unique的reads
# 过滤方法：MAPQ > 1
Rscript /lustre/user/liclab/lirf/Project/hic/bgi/5534N/code/creat.unique.r
sh ./tmp.code/uniquebam.sh

#此脚本用于去除PCR产生的duplication
#使用picard中的MarkDuplicates来做。
Rscript /lustre/user/liclab/lirf/Project/hic/bgi/5534N/code/creat.rmp.r
sh ./tmp.code/rmpbam.sh

#merge bamfiles& filted by chr(hap...Un_g) restriction loci
# 此处用于将2组分开比对的数据merge到一块
Rscript /lustre/user/liclab/lirf/Project/hic/bgi/5534N/code/creat.merge.r
sh ./tmp.code/merge.sh
####################################################################
# 到此为止，已经将原始的fastq文件经过比对产生SAM文件，转化为BAM文件，
# 按照染色体的位置排序，去除没有比对的序列， 去除PCR产生的duplicates
# 最后Merge到一块儿
####################################################################

# 此脚本根据reads的名字排序
Rscript /lustre/user/liclab/lirf/Project/hic/bgi/5534N/code/creat.sortname.r
sh ./tmp.code/sort.sh 

# 此脚本将BAM文件转化成SAM文件
Rscript /lustre/user/liclab/lirf/Project/hic/bgi/5534N/code/creat.bamtosam.r
sh ./tmp.code/bamtosam.sh

# 只提取SAM文件的前4列信息：名称，染色体信息，起始为止，终止位置
less ./tmp.data/BAM/*fq.sam | awk '{print $1,$2,$3,$4}' > ./tmp.data/BAM/sam_merge.txt

#此步过滤去除不在24条染色体上的reads: 得到文件 filter.txt
Rscript /lustre/user/liclab/lirf/Project/hic/bgi/5534N/code/creat.filter.py.r
sh ./tmp.code/filter.py.sh

# 此步过滤去 500bp以内无h3酶切位点的reads，得到文件 filter.hindIII.txt
Rscript /lustre/user/liclab/lirf/Project/hic/bgi/5534N/code/creat.hindIIIfilter500.py.r
sh ./tmp.code/hindIIIfilter500.py.sh

#此步将reads过滤完毕后Merge成reads pair在一行的文件， 得到文件 filter.hindIII.merge.txt
Rscript /lustre/user/liclab/lirf/Project/hic/bgi/5534N/code/creat.merge.py.r
sh ./tmp.code/merge.py.sh
################################################################
# 以上做完了reads 水平的过滤
################################################################

################################################################
#fragment水平的过滤
################################################################

# 此步只挑出位于同一条染色体的reads：
awk '$3 == $6{print}' ./tmp.data/BAM/filter.hindIII.merge.txt > ./tmp.data/BAM/cis.txt

# 此步对每对reads增加两列数据（1.[inward,outward,samestrand],2.gap size(fragement 之间的距离)）
python /lustre/user/liclab/lirf/Project/hic/bgi/5534N/code/cis.frage.py -i ./tmp.data/BAM/cis.txt -o ./tmp.data/BAM/cis.add.frage.txt

#去掉在同一个fragments上的reads.
less ./tmp.data/BAM/cis.add.frage.txt |awk '{if($9>=0){print}}' > ./tmp.data/BAM/cis.nosamfrag.txt

#画图Cutoff distance between read pairs 
mkdir Graph
mkdir ./Graph/heatmapForRawMarix
mkdir ./Graph/heatmapForNormalizedMarix
python /lustre/user/liclab/lirf/Project/hic/bgi/5534N/code/graph.py
python /lustre/user/liclab/lirf/Project/hic/bgi/5534N/code/sum.py
Rscript /lustre/user/liclab/lirf/Project/hic/bgi/5534N/code/plot.reads.R

#画图Error structure by distance
less ./tmp.data/BAM/cis.nosamfrag.txt |awk '{if($8==1 && $9>=0){print $9}}' > ./tmp.data/BAM/inward.txt
less ./tmp.data/BAM/cis.nosamfrag.txt |awk '{if($8==2 && $9>=0){print $9}}' > ./tmp.data/BAM/outward.txt
less ./tmp.data/BAM/cis.nosamfrag.txt |awk '{if($8==3 && $9>=0){print $9}}' > ./tmp.data/BAM/samestrand.txt
Rscript /lustre/user/liclab/lirf/Project/hic/bgi/5534N/code/plot.R 

# 此步做fragment过滤（inward>1kb,outword>25kb）：得到最后过滤完毕的reads
python /lustre/user/liclab/lirf/Project/hic/bgi/5534N/code/fragement_level_filter.py -i ./tmp.data/BAM/cis.nosamfrag.txt -o ./tmp.data/BAM/cis.filter.txt


###############################################################
#产生染色体内部相互作用矩阵并保存在./matrix/cis
###############################################################

mkdir ./tmp.data/result
mkdir ./matrix
mkdir ./matrix/cis
mkdir ./matrix/trans
mkdir ./matrix/cis/raw
mkdir ./matrix/trans/raw
mkdir ./matrix/cis/normalization
mkdir ./matrix/trans/normalization
for((i=1;i<=24;i++));
do
	less ./tmp.data/BAM/cis.filter.txt | awk "\$3 == $i && \$6 == $i {print}"  > ./tmp.data/result/chr${i}_reads.txt
#creat interaction matrix, 构建相互作用矩阵
	python /lustre/user/liclab/lirf/Project/hic/bgi/5534N/code/creat_cis_matrix.py -i ./tmp.data/result/chr${i}_reads.txt  -o ./matrix/cis/raw/chr${i}_200kb_rawmatrix.txt -chr $i -resolution 200000
done
###############################################################
#产生染色体之间的相互作用矩阵并保存在./matrix/trans
###############################################################

awk '$3 != $6{print}' ./tmp.data/BAM/filter.hindIII.merge.txt > ./tmp.data/BAM/trans.txt
for((i=1;i<=23;i++));
do
for((j=i+1;j<=24;j++));
do
  less ./tmp.data/BAM/trans.txt | awk "\$3 == $i && \$6 == $j||\$3 == $j && \$6 == $i {print}"  > ./tmp.data/result/chr${i}_${j}_reads.txt
#creat interaction matrix, 构建相互作用矩阵
	python /lustre/user/liclab/lirf/Project/hic/bgi/5534N/code/creat_trans_matrix.py -i ./tmp.data/result/chr${i}_${j}_reads.txt  -o ./matrix/trans/raw/chr${i}_${j}_200kb_rawmatrix.txt -chr1 $i -chr2 $j -resolution 200000
done
done

################################################################
#normalization 绘制heatmap
################################################################
# Rscript ./code/rcode_glm_normalization_cis.R &
# Rscript ./code/rcode_glm_normalization_trans.r &

#######################产生500kb 的相互作用矩阵并校正
sh /lustre/user/liclab/lirf/Project/hic/bgi/5534N/code/run_500kb.sh
cd ./resolution_500k/cis
Rscript /lustre/user/liclab/lirf/Project/hic/bgi/5534N/code/ice_500k.R

#######################产生40kb 的相互作用矩阵并校正

cd ../../
sh /lustre/user/liclab/lirf/Project/hic/bgi/5534N/code/run_40k.sh
cd ./resolution_40k/cis
Rscript /lustre/user/liclab/lirf/Project/hic/bgi/5534N/code/ice_40k.R

#######################Call TAD
Rscript /lustre/user/liclab/lirf/Project/hic/bgi/5534N/code/02.prepare_files.R
cd ./TAD_files
gzip *
cd ../
sh /lustre/user/liclab/lirf/Project/hic/bgi/5534N/code/Call_TAD_boundary_40kb.sh

