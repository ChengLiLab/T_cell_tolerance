mkdir resolution_10k
mkdir ./resolution_10k/cis
mkdir ./resolution_10k/trans
mkdir ./resolution_10k/cis/raw
mkdir ./resolution_10k/cis/normalization
mkdir ./resolution_10k/trans/raw
mkdir ./resolution_10k/trans/normalization


for((i=1;i<=24;i++));
do
    #less ./tmp.data/BAM/cis.filter.txt | awk "\$3 == $i && \$6 == $i {print}"  > ./tmp.data/result/chr${i}_reads.txt
    #creat interaction matrix, 构建相互作用矩阵
    python /lustre/user/liclab/lirf/Project/hic/bgi/5534N/code/creat_cis_matrix.py -i ./tmp.data/result/chr${i}_reads.txt  -o ./resolution_10k/cis/raw/chr${i}_10k_rawmatrix.txt -chr $i -resolution 10000
done
###############################################################
#产生染色体之间的相互作用矩阵并保存在./matrix/trans
###############################################################

#awk '$3 != $6{print}' ./tmp.data/BAM/filter.hindIII.merge.txt > ./tmp.data/BAM/trans.txt
# for((i=1;i<=23;i++));
# do
#     for((j=i+1;j<=24;j++));
#     do
#         #  less ./tmp.data/BAM/trans.txt | awk "\$3 == $i && \$6 == $j||\$3 == $j && \$6 == $i {print}"  > ./tmp.data/result/chr${i}_${j}_reads.txt
#         #creat interaction matrix, 构建相互作用矩阵
#         python /lustre/user/liclab/lirf/Project/hic/bgi/5534N/code/creat_trans_matrix.py -i ./tmp.data/result/chr${i}_${j}_reads.txt  -o ./resolution_20k/trans/raw/chr${i}_${j}_20k_rawmatrix.txt -chr1 $i -chr2 $j -resolution 20000
#     done
# done


