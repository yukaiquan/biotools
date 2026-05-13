#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os

input_bim = sys.argv[1]
output_bim = sys.argv[2]
output_map = sys.argv[3]

output_bim_open = open(output_bim, 'w')
res_list = []
res_dict = {}
with open(input_bim, 'r') as f:
    num = 0
    for line in f:
        if line.startswith('#'):
            continue
        line = line.strip()
        lines = line.split()
        chr_name = lines[0]
        if chr_name not in res_dict:
            num += 1
            res_dict[chr_name] = num
            output_bim_open.write(str(num) + '\t' + '\t'.join(lines[1:]) + '\n')
        else:
            output_bim_open.write(str(res_dict[chr_name]) + '\t' + '\t'.join(lines[1:]) + '\n')
            continue

output_bim_open.close()
with open(output_map, 'w') as f:
    for key, value in res_dict.items():
        f.write(key + '\t' + str(value) + '\n')
        
# python plinkreplacebim.py input.bim output.bim output.map