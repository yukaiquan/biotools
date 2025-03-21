#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# gatk must be installed and added to the PATH, and the gatk database must be created(version 4.5.0.0)
# tmp dir must be set to a large enough disk space and must be created before running the script
import sys
import os
import json

input_gatk_db = sys.argv[1]
# input_gatk_db = "GenomeNGSGATK"
input_chr_name = sys.argv[2]
# input_chr_name = "chr1A_1"
reference_fasta = sys.argv[3]
# reference_fasta = "chr1A_n.fasta"
tmp_dir = "tmp"
# output_file = "gatk_gvcf_split_output.tsv"
output_file = sys.argv[4]
input_window_size = 10000000
# input_step_size = 9000000
input_step_size = 10000000


# 检查文件是否存在
if not os.path.exists(tmp_dir):
    os.makedirs(tmp_dir)
if not os.path.exists(input_gatk_db):
    print("GATK database does not exist")
    sys.exit(1)
if not os.path.exists(reference_fasta):
    print("Reference fasta does not exist")
    sys.exit(1)

def read_db_chr(input_gatk_db, input_chr_name)->dict:
    gatk_db_vidmap_json = os.path.join(input_gatk_db, "vidmap.json")
    chr_len_dict = {}

    with open(gatk_db_vidmap_json, "r") as f:
        vidmap = json.load(f)
        for vid in vidmap:
            chr_list = vidmap[vid]
            for chr in chr_list:
                # print(chr)
                if chr['name'] == input_chr_name:
                    chr_len_dict[chr['name']] = chr['length']
    return chr_len_dict

# read db chr to a dict
chr_len_dict = read_db_chr(input_gatk_db, input_chr_name)
# print(chr_len_dict)

# generate bin list by chr length  eg:[0-10000000,9900000-19900000, ...]
def generate_bin_list(chr_len_dict,input_chr_name, input_window_size, input_step_size)->list:
    bin_list = []
    last_bin_end = 1
    length = int(chr_len_dict[input_chr_name])
    while last_bin_end < length and last_bin_end + input_window_size <= length:
        bin_list.append([last_bin_end, last_bin_end + input_window_size])
        last_bin_end += input_step_size
    if last_bin_end + input_window_size >= length:
        bin_list[-1][1] = length
    return bin_list

bin_list = generate_bin_list(chr_len_dict,input_chr_name, input_window_size, input_step_size)
# print(bin_list)

# gatk GenotypeGVCFs -OVI False -R chr1A_n.fasta --tmp-dir tmp -O chr1A_1.vcf.gz -V gendb://$my_database -L chr1A_1
# shell_list = []
shell_dict = {}
file_dict = {}

for bin in bin_list:
    file_name = f"{input_chr_name}_{bin[0]}_{bin[1]}.vcf.gz"
    shell = f"gatk GenotypeGVCFs -OVI False -R {reference_fasta} --tmp-dir {tmp_dir} -O {file_name} -V gendb://{input_gatk_db} -L {input_chr_name}:{bin[0]}-{bin[1]}"
    if file_name not in shell_dict:
        shell_dict[file_name] = shell
    else:
        print(f"file name {file_name} already exists")
        exit(1)
    if file_name not in file_dict:
        file_dict[file_name] = [input_chr_name, bin[0], bin[1]]
    else:
        print(f"file name {file_name} already exists")
        exit(1)

with open(output_file, "w") as f:
    for file_name in file_dict:
        vaules = file_dict[file_name]
        f.write(f"{file_name}\t{vaules[0]}\t{vaules[1]}\t{vaules[2]}\n")

for file_name in shell_dict.keys():
    value = shell_dict[file_name]
    print(f"{file_name}\t{value}")
            