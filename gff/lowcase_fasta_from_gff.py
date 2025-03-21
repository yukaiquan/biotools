#!/usr/bin/env python

import re
import sys

# 读取GFF文件并提取区域信息
def parse_gff(gff_file):
    regions = []
    with open(gff_file, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 5:
                chrom = parts[0]
                start = int(parts[3]) - 1  # GFF是1-based，Python是0-based
                end = int(parts[4])
                regions.append((chrom, start, end))
    return regions

# 将指定区域的碱基转换为小写
def convert_to_lowercase(fasta_sequences, regions):


    for chrom, start, end in regions:
        # 提取指定区域的序列
        region_chr = fasta_sequences[chrom]
        # 只更改指定区域，保持其他区域不变
        region_chr = region_chr[:start] + region_chr[start:end].lower() + region_chr[end:]
        # 将修改后的序列存储到字典中
        fasta_sequences[chrom] = region_chr

    return fasta_sequences

# 将修改后的序列写回到新的FASTA文件中
def write_to_new_fasta(modified_sequences, output_file):
    with open(output_file, 'w') as file:
        for chrom in modified_sequences:
            file.write(f">{chrom}\n")
            # file.write(f"{modified_sequences[chrom]}\n")
            # 每60个字符换行，保持FASTA格式
            for i in range(0, len(modified_sequences[chrom]), 60):
                file.write(f"{modified_sequences[chrom][i:i+60]}\n")

# 主程序
def main(gff_file, fasta_file, output_file):
    # 解析GFF文件
    regions = parse_gff(gff_file)
    
    # 假设FASTA文件的序列数据是按照染色体顺序存储的
    # 这里需要根据实际情况来调整，比如使用字典来存储每个染色体的序列
    # 以下代码仅为示例，需要根据实际情况调整
    fasta_sequences = {}  # 存储FASTA文件中的序列数据
    with open(fasta_file, 'r') as file:
        for line in file:
            if line.startswith('>'):
                chrom = line.strip()[1:]  # 假设FASTA头信息格式为">chromosome_name"
                fasta_sequences[chrom] = ''
            else:
                fasta_sequences[chrom] += line.strip()

    # 将指定区域的碱基转换为小写
    modified_sequences = convert_to_lowercase(fasta_sequences, regions)

    # 将修改后的序列写回到新的FASTA文件中
    write_to_new_fasta(modified_sequences, output_file)

if __name__ == "__main__":
    # gff_file = 'example.gff'  # GFF文件路径
    gff_file = sys.argv[1]  # GFF文件路径
    # fasta_file = 'example.fasta'  # FASTA文件路径
    fasta_file = sys.argv[2]  # FASTA文件路径
    # output_file = 'modified.fasta'  # 输出文件路径
    output_file = sys.argv[3]  # 输出文件路径
    main(gff_file, fasta_file, output_file)