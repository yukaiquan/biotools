#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import gzip


input = sys.argv[1]
output = sys.argv[2]

header = ""
res_dict = {}

# with open(input, "r") as f:
with gzip.open(input, "rt") as f:
    # for line in f:
    for line in f:
        if line.startswith("#"):
            header += line
        else:
            lines = line.strip().split("\t")
            chr_name = lines[0]
            if chr_name == "cpDNA":
                continue
            if chr_name not in res_dict:
                res_dict[chr_name] = [line]
            else:
                res_dict[chr_name].append(line)

# 除使用位置外使用优先级gene mRNA exon CDS的顺序来排序
# sorted_tem = {"gene": 0, "mRNA": 1, "exon": 2, "CDS": 3}
sorted_tem_strand_1 = {"gene": 0, "mRNA": 1, "five_prime_UTR":2,"exon": 3, "CDS": 4,"three_prime_UTR":5}
sorted_tem_strand_2 = {"gene": 0, "mRNA": 1, "three_prime_UTR":2,"exon": 3, "CDS": 4, "five_prime_UTR":5}


for chr_name in res_dict:
    content = res_dict[chr_name]
    #首先根据坐标排序,再根据优先级排序
    # res_dict[chr_name] = sorted(content, key=lambda x: (int(x.split("\t")[3]), sorted_tem_strand_1[x.split("\t")[2]]))
    # 更智能的排序方式，正链用sorted_tem_strand_1，负链用sorted_tem_strand_2
    res_dict[chr_name] = sorted(content, key=lambda x: (int(x.split("\t")[3]), sorted_tem_strand_1[x.split("\t")[2]]) if x.split("\t")[6] == "+" else (int(x.split("\t")[3]), sorted_tem_strand_2[x.split("\t")[2]]))

# 按照染色体顺序排序
chr_list = sorted(res_dict.keys())

with open(output, "w") as f:
    f.write(header)
    for chr_name in chr_list:
        for line in res_dict[chr_name]:
            f.write(line)
            
