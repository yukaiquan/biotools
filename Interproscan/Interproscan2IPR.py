#!/usr/bin/env python
# coding: utf-8
import sys
import os
import argparse
import subprocess


def arg_parse():
    parser = argparse.ArgumentParser(description="Parse Interproscan Results to IPR Table")
    parser.add_argument("-i", "--input", help="Input Interproscan Results", required=True)
    parser.add_argument("-o", "--output", help="Output IPR Table", required=True)
    args = parser.parse_args()
    return args

# # input_ipr_tsv = sys.argv[1]
# input_ipr_tsv = "./test_ipr.tsv"
# # output_ipr_tsv = sys.argv[2]
# output_ipr_tsv = "./test_ipr_out.tsv"

args = arg_parse()
input_ipr_tsv = args.input
output_ipr_tsv = args.output


outpu_ann_dict = {}
with open(input_ipr_tsv, 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        else:
            lines = line.split('\t')
            gene_id = lines[0]
            pep_len = lines[2]
            align_start = lines[6]
            align_end = lines[7]
            e_value = lines[8]
            ipr_id = lines[11]
            ipr_desc = lines[12]
            if "IPR" in ipr_id:
                if gene_id not in outpu_ann_dict:
                    outpu_ann_dict[gene_id] = [(ipr_id, ipr_desc, align_start, align_end, pep_len, e_value)]
                    
                else:
                    outpu_ann_dict[gene_id].append((ipr_id, ipr_desc, align_start, align_end, pep_len, e_value))

with open(output_ipr_tsv, 'w') as f:
    f.write("#GeneID\tIPR_ID\tIPR_DESC\tAlign_start\tAlign_end\tPeptide_length\tE_value\n")
    for gene_id in outpu_ann_dict:
        ipr_ann = outpu_ann_dict[gene_id]
        for ipr_info in ipr_ann:
            f.write(gene_id + "\t" + "\t".join(ipr_info) + "\n")