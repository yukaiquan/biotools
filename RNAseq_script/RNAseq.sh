#!/bin/bash 
#SBATCH -J salt_rnaseq
#SBATCH -p hebhcnormal01
#SBATCH -n 32
#SBATCH -o std.out.%j
#SBATCH -e std.err.%j

# 清除当前环境加载的软件
module purge

# 激活当前环境中的环境变量
source /public/home/acfaa2ssz7/.bashrc
conda activate bwa

python RNAseq.py -s sample.txt -r /public/home/acfaa2ssz7/genome/SFS_STAR/ -i /public/home/acfaa2ssz7/genome/bwaindex/sfs/ACD.sativa.fasta -g /public/home/acfaa2ssz7/genome/ACD_sativa_sfs.gtf -t 32
wait
