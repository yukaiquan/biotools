#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os


input = sys.argv[1]
output = sys.argv[2]

header = ""
res_dict = {}

with open(input, "r") as f:
    for line in f:
        if line.startswith("#"):
            header += line
        else:
            lines = line.strip().split("\t")
            chr_name = lines[0]
            if chr_name not in res_dict:
                res_dict[chr_name] = [line]
            else:
                res_dict[chr_name].append(line)

# 除使用位置外使用优先级gene mRNA exon CDS的顺序来排序
sorted_tem = {"gene": 0, "mRNA": 1, "exon": 2, "CDS": 3}
for chr_name in res_dict:
    #首先根据坐标排序,再根据优先级排序
    res_dict[chr_name] = sorted(res_dict[chr_name], key=lambda x: (int(x.split("\t")[3]), sorted_tem[x.split("\t")[2]]))

# 按照染色体顺序排序
chr_list = sorted(res_dict.keys())

with open(output, "w") as f:
    f.write(header)
    for chr_name in chr_list:
        for line in res_dict[chr_name]:
            f.write(line)
            