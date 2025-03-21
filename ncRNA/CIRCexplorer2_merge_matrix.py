#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys



input_cir2_list_file = sys.argv[1]
output_file_prefix = sys.argv[2]
output_folder = sys.argv[3]

file_list = []
with open(input_cir2_list_file, 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        file_list.append(line.strip())


count_dict = {}

total_count = 0

for file_name in file_list:
    with open(file_name, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            line_list = line.strip().split('\t')
            chr_name = line_list[0]
            start = line_list[1]
            end = line_list[2]
            key = chr_name + '_' + start + '_' + end
            count = line_list[12]
            total_count += int(count)
            if key in count_dict:
                count_dict[key].append((file_name, count))
            else:
                count_dict[key] = [(file_name, count)]

count_matrix = []
srpbm_matrix = []

for key in count_dict.keys():
    count_list = count_dict[key]
    row_list = [key]
    row_srpbm_list = [key]
    poss = key.split('_')
    start = poss[1]
    end = poss[2]
    gene_len = int(end) - int(start) + 1
    for file_name in file_list:
        if file_name in [x[0] for x in count_list]:
            value = 0
            for count_tuple in count_list:
                if count_tuple[0] == file_name:
                    value = int(count_tuple[1])
                    break
            srpbm_value = (float(value) * 1000000000) / (total_count * gene_len)
            row_list.append(str(value))
            row_srpbm_list.append(str(srpbm_value))
        else:
            row_list.append('0')
            row_srpbm_list.append('0')
            
    count_matrix.append(row_list)
    srpbm_matrix.append(row_srpbm_list)

out_put_count_file = output_folder + '/' + output_file_prefix + '_count.txt'
out_put_srpbm_file = output_folder + '/' + output_file_prefix + '_srpbm.txt'

head_list = [x.split('/')[-1] for x in file_list]
with open(out_put_count_file, 'w') as f:
    f.write('ID' + '\t' + '\t'.join(head_list) + '\n')
    for i in range(len(count_matrix)):
        f.write('\t'.join(count_matrix[i]) + '\n')

with open(out_put_srpbm_file, 'w') as f:
    f.write('ID' + '\t' + '\t'.join(head_list) + '\n')
    for i in range(len(srpbm_matrix)):
        f.write('\t'.join(srpbm_matrix[i]) + '\n')



