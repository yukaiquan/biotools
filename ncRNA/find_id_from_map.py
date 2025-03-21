#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys


input_map_file = sys.argv[1]
input_file_list = sys.argv[2]
sep_file = sys.argv[3]
output_map_file = sys.argv[4]

map_dict = {}

with open(input_map_file, 'r') as f:
    for line in f:
        line = line.strip()
        if line:
            key, value = line.split(sep_file)
            if key in map_dict:
                map_dict[key].append(value)
            else:
                map_dict[key] = [value]

output_map_file_open = open(output_map_file, 'w')

with open(input_file_list, 'r') as f:
    for line in f:
        line = line.strip()
        if line:
            key = line
            if key in map_dict:
                for value in map_dict[key]:
                    output_map_file_open.write(key + sep_file + value + '\n')

output_map_file_open.close()

