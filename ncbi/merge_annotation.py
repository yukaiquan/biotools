#!/usr/bin/env python3
import sys
import os


input_nr = sys.argv[1]
input_kog = sys.argv[2]
input_swiss_prot = sys.argv[3]
output_file = sys.argv[4]


result_id_list: list = []


nr_ann_dict: dict = {}
with open(input_nr, "r") as f:
    for line in f:
        line = line.strip().split("\t")
        gene_id = line[0]
        nr_id = line[1]
        result_id_list.append(gene_id)
        if nr_id in nr_ann_dict:
            nr_ann_dict[gene_id].append(nr_id)
        else:
            nr_ann_dict[gene_id] = [nr_id]

kog_ann_dict: dict = {}
with open(input_kog, "r") as f:
    for line in f:
        line = line.strip().split("\t")
        gene_id = line[0]
        kog_id = line[2]
        result_id_list.append(gene_id)
        if kog_id in kog_ann_dict:
            kog_ann_dict[gene_id].append(kog_id)
        else:
            kog_ann_dict[gene_id] = [kog_id]

swiss_prot_ann_dict: dict = {}
with open(input_swiss_prot, "r") as f:
    for line in f:
        line = line.strip().split("\t")
        gene_id = line[0]
        swiss_prot_id = line[4]
        result_id_list.append(gene_id)
        if swiss_prot_id in swiss_prot_ann_dict:
            swiss_prot_ann_dict[gene_id].append(swiss_prot_id)
        else:
            swiss_prot_ann_dict[gene_id] = [swiss_prot_id]

# 对result_id_list去重排序
result_id_list = list(set(result_id_list))
result_id_list.sort()

# 输出文件
with open(output_file, "w") as f:
    for gene_id in result_id_list:
        f.write(gene_id + "\t")
        if gene_id in nr_ann_dict:
            f.write(";".join(nr_ann_dict[gene_id]) + "\t")
        else:
            f.write("NA\t")
        if gene_id in kog_ann_dict:
            f.write(";".join(kog_ann_dict[gene_id]) + "\t")
        else:
            f.write("NA\t")
        if gene_id in swiss_prot_ann_dict:
            f.write(";".join(swiss_prot_ann_dict[gene_id]) + "\t")
        else:
            f.write("NA\t")
        f.write("\n")
