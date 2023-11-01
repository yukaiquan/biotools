#!/usr/bin/env python3

import os
import sys

file_input = "CDvsACD.matrix"

marix_list:list = []

with open(file_input, "r") as f:
    for line in f:
        line = line.strip()
        # 分组排序
        lines = line.split("\t")
        a_name = lines[0]
        a_start = int(lines[1])
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
        marix_list_new.append(marix_list[i])
    else:
        if marix_list[i][0] == marix_list[i - 1][0] and marix_list[i][1] - marix_list[i - 1][2] > length:
            if marix_list[i][0] == marix_list[i - 1][0] and marix_list[i+1][1] - marix_list[i][2] > length:
                continue
            else:
                marix_list_new.append(marix_list[i])
        else:
            marix_list_new.append(marix_list[i])
    
