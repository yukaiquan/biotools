#!/usr/bin/env python3

import os
import sys

interproscan_input = sys.argv[1]
interproscan_output = sys.argv[2]

ipr = {}
gene_ipr = {}

with open(interproscan_input, 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        line_split = line.split('\t')
        gene = line_split[0]
        ipr_id = line_split[11]
        ipr_desc = line_split[12]
        if ipr_id.startswith('IPR'):
            if gene not in gene_ipr:
                gene_ipr[gene] = [ipr_id]
            else:
                gene_ipr[gene].append(ipr_id)
            if ipr_id not in ipr:
                ipr[ipr_id] = ipr_desc
            else:
                continue

with open(interproscan_output, 'w') as f:
    f.write('Gene\tIPR\tDescription\n')
    for gene in gene_ipr.keys():
        value = gene_ipr[gene]
        # 去重
        value = list(set(value))
        for ipr_id in value:
            ipr_desc = ipr[ipr_id]
            f.write(f'{gene}\t{ipr_id}\t{ipr_desc}\n')
