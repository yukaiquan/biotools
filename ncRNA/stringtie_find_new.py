#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os

input_stringtie_gff = sys.argv[1]
output_stringtie_gff = sys.argv[2]

new_gene = False

output_open = open(output_stringtie_gff, 'w')
new_gff_list = []
new_gff_dict = {}

class_code = ""
length_tr = 0

with open(input_stringtie_gff, 'r') as f:
    for line in f:
        if line.startswith("#"):
            continue
        else:
            line_split = line.strip().split("\t")
            if line_split[2] == "transcript":
                if 'class_code' in line_split[8]:
                    # 正则提取class_code=u
                    class_code = line_split[8].split("class_code")[1].split(";")[0].strip().replace('"', '')
                start = line_split[3]
                end = line_split[4]
                length_tr = int(end) - int(start) + 1
            # if class_code == "u" or class_code == "x" or class_code == "j" or class_code == "i" or class_code == "o":
            if class_code == "u" or class_code == "x":
                # output_open.write(line)
                # if length_tr >= 200:
                    tr_name = line_split[8].split(";")[0].replace("transcript_id ", "").strip().replace('"', '')
                    if tr_name not in new_gff_dict:
                        new_gff_dict[tr_name] = [line]
                    else:
                        new_gff_dict[tr_name].append(line)
            else:
                continue




for key in new_gff_dict.keys():
    value = new_gff_dict[key]
    if len(value) >2:
        for i in value:
            output_open.write(i)
                

output_open.close()
