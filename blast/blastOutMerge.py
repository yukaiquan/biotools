#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import os


# This script merges the blast output files into one file for each query
# Usage: python blastOutMerge.py <blastOutTfile> <blastOutTfile>

# Read in the blast output files
blastOutTfile = sys.argv[1]
blastOutQfile = sys.argv[2]

# 读取blast输出文件
blastOut1 = {}
blastOut2 = {}

# 读取blast输出文件
with open(blastOutTfile, 'r') as f:
    for line in f:
        line = line.strip()
        if line == '':
            continue
        line = line.split('\t')
        if line[0] not in blastOut1:
            blastOut1[line[0]] = []
        blastOut1[line[0]].append(line[1])

with open(blastOutQfile, 'r') as f:
    for line in f:
        line = line.strip()
        if line == '':
            continue
        line = line.split('\t')
        if line[1] not in blastOut2:
            blastOut2[line[1]] = []
        blastOut2[line[1]].append(line[0])

# Merge the blast output files
blastOut = {}
for key in blastOut1:
    if key not in blastOut:
        blastOut[key] = []
    blastOut[key] = blastOut1[key]

# 输出文件
with open('blastOutMerge.txt', 'w') as f:
    for key in blastOut:
        f.write(key + '\t' + ','.join(blastOut[key]) + '\n')
