#!/bin/bash 
#SBATCH -J rnaseq
#SBATCH -p xahdtest
#SBATCH -n 32
#SBATCH --gres=dcu:4	#指定每个节点使用4块DCU卡
#SBATCH -o std.out.%j
#SBATCH -e std.err.%j

# 清除当前环境加载的软件
module purge

# 激活当前环境中的环境变量
source /work/home/acfaa2ssz7/.bashrc
conda activate bwa

python ./RNAseq.py sample.txt 01_fastp 02_bam 03_feature 04_matrix 05_add 06_dup 07_gvcf 08_vcf /work/home/acfaa2ssz7/genome/bwaindex/sfs/ACD.sativa.fasta /work/home/acfaa2ssz7/genome/star/sfs/SFS_STAR /work/home/acfaa2ssz7/genome/star/sfs/ACD_sativa_sfs.gtf 16
wait
