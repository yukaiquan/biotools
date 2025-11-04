#!/usr/bin/env python

import sys
# 从wgdi生成的gff提取最长的转录本代表该基因

input_gff = sys.argv[1]
input_pep = sys.argv[2]
output_pep = sys.argv[3]

gene_len = {}
seq_dict = {}

with open(input_gff, 'r') as gff_file:
    for line in gff_file:
        if line.startswith("#"):
            continue
        line_split = line.split("\t")
        gene_name = line_split[1]
        start = int(line_split[2])
        end = int(line_split[3])
        tr_name = line_split[6].strip()
        if gene_name not in gene_len:
            gene_len[gene_name] = [tr_name, start, end]
        else:
            if abs(end - start) > abs(gene_len[gene_name][2] - gene_len[gene_name][1]):
                gene_len[gene_name] = [tr_name, start, end]

with open(input_pep, 'r') as pep_file:
    for line in pep_file:
        if line.startswith(">"):
            tr_name = line.strip().split(" ")[0][1:].strip()
            # print(tr_name)
            seq_dict[tr_name] = ""
        else:
            seq_dict[tr_name] += line.strip()

with open(output_pep, 'w') as output_file:
    for gene_name in gene_len.keys():
        tr_name = gene_len[gene_name][0].strip()
        if tr_name not in seq_dict:
            print(tr_name + " not found in pep file")
            continue
        output_file.write(">" + gene_name + "\n")
        # output_file.write(seq_dict[tr_name])
        seq = seq_dict[tr_name].replace(".", "").replace("*", "")
        # 一行60个字符
        for i in range(0, len(seq), 60):
            output_file.write(seq[i:i+60] + "\n")
