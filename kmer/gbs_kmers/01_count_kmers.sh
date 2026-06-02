#!/bin/bash
# ============================================================================
# Script: 01_count_kmers.sh
# Description: Count k-mers for GBS (Genotyping-by-Sequencing) samples
# ============================================================================

# ----------------------------------------------------------------------------
# Step 1: Submit batch jobs on SLURM cluster to calculate k-mer counts for 666 GBS samples
# ----------------------------------------------------------------------------

# Define path to the SNP BBDuk script
SNPBBDUK="/publicssd/share/h13713/soft/bbmap/snp_bbduk.zsh"

# Generate batch submission file
# Iterate through all paired-end FASTQ files (*_2.f*q.gz)
# Check if output file already exists, if not, add to submission list
ls */*_2.f*q.gz | while read id; do
    # Extract sample name from file path
    name=${id%%/*}
    
    # Define expected output file path
    output_file="${name}/${name}_wgs_count.tsv.gz"
    
    # Only submit jobs for samples without existing output
    if [ ! -f "$output_file" ]; then
        # Format: sample_name \t command
        echo -e "${name}\tbash ${SNPBBDUK} GBS_uniq.fasta ${name}"
    fi
done | sort | uniq > sbatch.tsv

# Submit jobs to SLURM cluster
# Parameters:
#   -i: input file with sample names and commands
#   -p: partition name (com300)
#   -m: memory allocation (20G per job)
#   -N: number of nodes (1)
#   -n: number of CPUs per task (4)
#   -e: environment setup (bwa)
python /publicssd/share/h13713/soft/sbatch_script.py -i sbatch.tsv -p com300 -m 20G -N 1 -n 4 -e bwa
# https://github.com/yukaiquan/biotools/blob/main/slurm/sbatch_script.py

# ----------------------------------------------------------------------------
# Step 2: Extract k-mer count tables from original SNP matrix for 9000 GBS samples
# ----------------------------------------------------------------------------

# Convert GBS VCF to GBS dataset
# This step processes the original VCF file to generate k-mer counts
#
# Logic:
#   - Genotypes (GT) are used to distinguish homozygous (homo) vs heterozygous (hete)
#   - DP (total depth) is used for counting
#   - For homozygous samples (0/0 or 1/1): record DP directly
#   - For heterozygous samples (0/1): record both DP and DV (variant depth)
#
# Output:
#   - One folder per sample
#   - Each folder contains: {sample}_gbs_count.tsv.gz

python getkmerbyvcf.py GBS_uniq.fasta all.gz


