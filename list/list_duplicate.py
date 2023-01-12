#!/use/bin/env python3
# -*- coding: utf-8 -*-
# Author: yukaiquan
# Email: 1962568272@qq.com
# usage: python3 list_duplicate.py input_list output_file colums
# colums is duplicate, split is tab
import sys

input_list: str = sys.argv[1]
output_file: str = sys.argv[2]
colmes: int = int(sys.argv[3])

list_file: list = []
dup_list: list = []

with open(input_list, 'r') as f:
    line = f.readlines()
    for l in line:
        l = l.strip().split('\t')
        if l[colmes] not in dup_list:
            dup_list.append(l[colmes])
            list_file.append(l)
        else:
            print(l)

with open(output_file, 'w') as f:
    for l in list_file:
        f.write('\t'.join(l) + '\n')
