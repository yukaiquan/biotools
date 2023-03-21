#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse as ap

# input_status = 'BY2_MERGEBAM_STAT.txt'
# output_table = 'status_table.txt'
# 读取参数
parser = ap.ArgumentParser(description='bam status to table')
parser.add_argument(
    '-i', '--input', help='input bam status file', required=True)
parser.add_argument('-o', '--output', help='output table file', required=True)
args = parser.parse_args()
input_status = args.input
output_table = args.output

# SFS Genome_length
Genome_length = 10757474745

out_list = []
with open(input_status, 'r') as f:
    lines = f.readlines()
    for line in lines:
        if line.startswith('#'):
            continue
        else:
            if line.startswith('SN'):
                line = line.strip().split('\t')
                # depth=bases mapped/genome length
                if line[1] == 'bases mapped:':
                    bases_mapped = line[2]
                    depth = float(bases_mapped) / Genome_length
                    out_list.append(['bases mapped', bases_mapped])
                    out_list.append(['depth', depth])
                else:
                    sample = line[1].replace(':', '').strip()
                    status = line[2].strip()
                    out_list.append([sample, status])
            else:
                continue

with open(output_table, 'w') as f:
    for i in out_list:
        f.write('\t'.join(map(str, i)) + '\n')
