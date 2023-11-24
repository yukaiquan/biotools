#!/usr/bin/env python3
import numpy as np
import os
import sys


# file_input = "CdvsAcd.txt"
file_input = sys.argv[1]
# chr_length = "Chr_length.txt"
chr_length = sys.argv[2]
chr_status = "CdvsAcd_status.txt"
block_re = "CdvsAcd_block.txt"
SIZE = 1000000

# 利用numpy生成区间列表


def generate_interval_list(genome_chr_len, window_size) -> list:
    # interval_list = np.arange(0, genome_chr_len + 1, window_size)
    # 最后一个区间的长度可能不足1Mb，所以需要单独处理
    interval_list = np.arange(0, genome_chr_len, window_size)
    interval_list = np.append(interval_list, genome_chr_len)
    interval_list = interval_list.astype(str)
    return interval_list


# 根据区间列表生成区间字典
def generate_interval_dict(chr_name: str, interval_list: list) -> dict:
    chr_dict = {}
    chr_dict[chr_name] = {}
    for i in range(len(interval_list) - 1):
        # interval_dict = {}
        # interval_dict[interval_list[i] + '-' + interval_list[i + 1]] = []
        # chr_dict[chr_name] = interval_dict
        # 嵌套字典
        chr_dict[chr_name][interval_list[i] + '-' + interval_list[i + 1]] = []

    return chr_dict


def parase_dict(length_dict: dict, window_size: int) -> dict:
    all_chr_dict = {}
    for key in length_dict.keys():
        chr_name = key
        chr_length = length_dict[key]
        interval_list = generate_interval_list(chr_length, window_size)
        chr_dict = generate_interval_dict(chr_name, interval_list)
        # 合并至大字典
        all_chr_dict.update(chr_dict)
    return all_chr_dict

# 根据数值解析出区间的位置


def caculate_interval_from_value(size: int, value: int, max_value: int) -> str:
    if value < size:
        return '0-' + str(size)
    else:
        # 除数取整
        div = value // size
        start_value = div * size
        end_value = start_value + size
        if end_value > max_value:
            end_value = max_value
        # print(str(start_value) + '-' + str(end_value))
        return str(start_value) + '-' + str(end_value)


marix_list: list = []

with open(file_input, "r") as f:
    for line in f:
        line = line.strip()
        # 分组排序
        lines = line.split("\t")
        a_name = lines[0]
        a_start = int(lines[1])
        # print(lines)
        a_end = int(lines[2])
        if a_start > a_end:
            a_start, a_end = a_end, a_start
        b_name = lines[3]
        b_start = int(lines[4])
        b_end = int(lines[5])
        if b_start > b_end:
            b_start, b_end = b_end, b_start
        marix_list.append([a_name, a_start, a_end, b_name, b_start, b_end])

# 根据a_name和a_start排序
marix_list.sort(key=lambda x: (x[0], x[1]))
# 间隔2MB内没有则去除
length = 2000000
marix_list_new = []
for i in range(len(marix_list)):
    if i == 0:
        # marix_list_new.append(marix_list[i])
        if marix_list[i+1][1] - marix_list[i][2] > length:
            continue
        else:
            marix_list_new.append(marix_list[i])
    else:
        if marix_list[i][0] == marix_list[i - 1][0] and marix_list[i][1] - marix_list[i - 1][2] > length:
            if i + 1 == len(marix_list):
                continue
            # print(marix_list[i])
            # print(i)
            if marix_list[i][0] == marix_list[i - 1][0] and marix_list[i+1][1] - marix_list[i][2] > length:
                continue
            else:
                marix_list_new.append(marix_list[i])
        else:
            if marix_list[i][0] != marix_list[i - 1][0] and marix_list[i][1] - marix_list[i-1][2] > length:
                continue
            else:
                marix_list_new.append(marix_list[i])

# 拆分block块，每1MB一个，统计每MB出现的reads数量
marix_list_new_block = []

A_chr_length: dict = {}
B_chr_length: dict = {}
with open(chr_length, 'r') as f:
    for line in f:
        lines = line.strip().split('\t')
        # print(lines)
        if lines[0] == 'A':
            if lines[1] not in A_chr_length:
                A_chr_length[lines[1]] = int(lines[2])
        elif lines[0] == 'B':
            if lines[1] not in B_chr_length:
                B_chr_length[lines[1]] = int(lines[2])

A_interval_dict = parase_dict(A_chr_length, SIZE)
B_interval_dict = parase_dict(B_chr_length, SIZE)

# 分配区间
for i in marix_list_new:
    chr_A_name = i[0]
    chr_A_start = i[1]
    chr_A_end = i[2]
    chr_B_name = i[3]
    chr_B_start = i[4]
    chr_B_end = i[5]
    if chr_A_name == 'chrUn' or chr_B_name == 'chrUn':
        continue
    A_interval = caculate_interval_from_value(
        SIZE, chr_A_start, A_chr_length[chr_A_name])
    B_interval = caculate_interval_from_value(
        SIZE, chr_B_start, B_chr_length[chr_B_name]
    )
    if chr_A_name in A_interval_dict:

        if A_interval in A_interval_dict[chr_A_name]:
            A_interval_dict[chr_A_name][A_interval].append(
                [B_interval, chr_A_name, chr_A_start, chr_A_end, chr_B_name, chr_B_start, chr_B_end])
        else:
            # print(chr_A_name, interval)
            # print(SIZE, chr_A_start, A_chr_length[chr_A_name])
            print('error', chr_A_name, A_interval)
            continue
    if chr_B_name in B_interval_dict:
        if B_interval in B_interval_dict[chr_B_name]:
            B_interval_dict[chr_B_name][B_interval].append(
                [A_interval, chr_A_name, chr_A_start, chr_A_end, chr_B_name, chr_B_start, chr_B_end])
        else:
            print('error', chr_B_name, B_interval)
            continue

# 统计染色体层面的基因共线性块数量
chr_block_num = {}

# 生成区间对应坐标关系
block_xy_num = {}

for chr_A_name in A_interval_dict:
    if chr_A_name == 'chrUn':
        continue
    if chr_A_name not in chr_block_num:
        chr_block_num[chr_A_name] = {}
    for interval in A_interval_dict[chr_A_name]:
        if len(A_interval_dict[chr_A_name][interval]) == 0:
            continue

        block_name = chr_A_name + '_' + interval
        if block_name not in block_xy_num:
            block_xy_num[block_name] = {}
        for value in A_interval_dict[chr_A_name][interval]:
            chr_B_name = value[4]
            if chr_B_name not in chr_block_num[chr_A_name]:
                chr_block_num[chr_A_name][chr_B_name] = 0
            else:
                chr_block_num[chr_A_name][chr_B_name] += 1
            block_name_b = chr_B_name + '_' + value[0]
            if block_name_b not in block_xy_num[block_name]:
                block_xy_num[block_name][block_name_b] = 0
            else:
                block_xy_num[block_name][block_name_b] += 1

with open(chr_status, 'w') as w:
    for i in chr_block_num:
        for j in chr_block_num[i]:
            w.write(i + '\t' + j + '\t' + str(chr_block_num[i][j]) + '\n')

with open(block_re, 'w') as w:
    for i in block_xy_num:
        for j in block_xy_num[i]:
            w.write("\t".join(i.split('_')) + '\t' + "\t".join(j.split('_')
                                                               ) + '\t' + str(block_xy_num[i][j]) + '\n')
