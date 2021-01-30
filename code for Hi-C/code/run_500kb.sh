mkdir resolution_500k
mkdir ./resolution_500k/cis
mkdir ./resolution_500k/cis/raw

for((i=1;i<=24;i++));
do
    #less ./tmp.data/BAM/cis.filter.txt | awk "\$3 == $i && \$6 == $i {print}"  > ./tmp.data/result/chr${i}_reads.txt
    #creat interaction matrix, 构建相互作用矩阵
    python /lustre/user/liclab/lirf/Project/hic/bgi/5534N/code/creat_cis_matrix.py -i ./tmp.data/result/chr${i}_reads.txt  -o ./resolution_500k/cis/raw/chr${i}_500k_rawmatrix.txt -chr $i -resolution 500000
done

mkdir ./resolution_500k/trans
mkdir ./resolution_500k/trans/raw



for((i=1;i<=23;i++));
do
    for((j=i+1;j<=24;j++));
    do
        #  less ./tmp.data/BAM/trans.txt | awk "\$3 == $i && \$6 == $j||\$3 == $j && \$6 == $i {print}"  > ./tmp.data/result/chr${i}_${j}_reads.txt
        #creat interaction matrix, 构建相互作用矩阵
        python /lustre/user/liclab/lirf/Project/hic/bgi/5534N/code/creat_trans_matrix.py -i ./tmp.data/result/chr${i}_${j}_reads.txt  -o ./resolution_500k/trans/raw/chr${i}_${j}_500k_rawmatrix.txt -chr1 $i -chr2 $j -resolution 500000
    done
done