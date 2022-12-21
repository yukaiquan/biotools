#!/bin/bash 
#SBATCH -J bwa_gatk
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

srun fastp -i /work/home/acfaa2ssz7/YKQNGS/sfs/01_fastq/V300030890_L04_2_1.fq.gz -I /work/home/acfaa2ssz7/YKQNGS/sfs/01_fastq/V300030890_L04_2_2.fq.gz \
		-o /work/home/acfaa2ssz7/YKQNGS/sfs/01_fastq/out_V300030890_L04_2_1.fq.gz -O /work/home/acfaa2ssz7/YKQNGS/sfs/01_fastq/out_V300030890_L04_2_2.fq.gz -w 32
