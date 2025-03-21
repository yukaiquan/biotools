#!/usr/bin/env python
# -*- coding: utf-8 -*-
from tqdm import tqdm
import os
import sys
import time


input_evm_mrna_include = sys.argv[1]
input_evm_all = sys.argv[2]
output_file = sys.argv[3]


includ_dict = []

with open(input_evm_mrna_include, 'r') as f:
    for line in tqdm(f, "Reading mRNA include file..."):
        if line.startswith('#'):
            continue
        else:
            lines = line.strip().split('\t')
            mrna_id = lines[8].split(';')[0].split('=')[1]
            gene_id = lines[8].split(';')[1].split('=')[1]
            includ_dict.append(mrna_id)
            includ_dict.append(gene_id)

out_list = []

with open(input_evm_all, 'r') as f:
        for line in tqdm(f, "Removing TEs..."):
            if line.startswith('#'):
                # out_open.write(line)
                out_list.append(line)
            else:
                 lines = line.strip().split('\t')
                 type_name = lines[2]
                 if type_name == 'mRNA':
                    mrna_id = lines[8].split(';')[0].split('=')[1]
                    if mrna_id in includ_dict:
                        # out_open.write(line)
                        out_list.append(line)
                    else:
                        continue
                 elif type_name == 'gene':
                    gene_id = lines[8].split(';')[0].split('=')[1]
                    if gene_id in includ_dict:
                        # out_open.write(line)
                        out_list.append(line)
                    else:
                        continue
                 else:
                     mrna_id = lines[8].split(';')[1].split('=')[1]
                     if mrna_id in includ_dict:
                        #  out_open.write(line)
                        out_list.append(line)
                     else:
                         continue
                    
                     

with open(output_file, 'w') as f:
    for line in tqdm(out_list, "Writing output file..."):
        f.write(line)
