#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys


input_gene_count_file = sys.argv[1]
input_gtf_file = sys.argv[2]
output_file = sys.argv[3]


gene_count_dict = {}
count = 3

with open(input_gene_count_file, "r") as f:
    for line in f:
        if line.startswith("transcript_id"):
            continue
        lines = line.strip().split(",")
        tr_name = lines[0]
        count_len = 0
        for i in range(1, len(lines)):
            if int(lines[i]) >= count:
                count_len += 1
        if count_len == len(lines) - 1:
            gene_count_dict[tr_name] = 1
        else:
            gene_count_dict[tr_name] = 0

out_open = open(output_file, "w")

with open(input_gtf_file, "r") as f:
    for line in f:
        if line.startswith("#"):
            out_open.write(line)
            continue
        lines = line.strip().split("\t")
        tr_name = lines[8].split(";")[0].split("transcript_id ")[1].strip().replace('"', "")
        if tr_name in gene_count_dict:
            if gene_count_dict[tr_name] == 1:
                out_open.write(line)
                
out_open.close()