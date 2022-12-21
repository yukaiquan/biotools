#!/bin/bash 
#SBATCH -J fastp
#SBATCH -p xahdtest
#SBATCH -n 32
#SBATCH -N 1
#SBATCH --ntasks-per-node=32
#SBATCH -o std.out.%j
#SBATCH -e std.err.%j

# 清除当前环境加载的软件
module purge

# 激活当前环境中的环境变量
source /work/home/acfaa2ssz7/.bashrc
conda activate bwa

python RNAseq.py sample.txt 01_fastp 16 2