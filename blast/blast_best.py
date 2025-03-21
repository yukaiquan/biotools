#!/usr/bin/python
# -*- coding: utf-8 -*-
# conda activate python3
"""
    作者：徐诗芬
    内容：读取并输出所有blast结果中最好的第一条hit结果
    日期：2021.1.22
"""
import sys

def usage():
    print('Usage: python blast_best.py [input_file] [outfile]')


def main():

    # global name
    inf = open(sys.argv[1], 'r')

    out_list = []
    flag_list = []
    for line in inf:
        line = line.strip()
        name = line.split("\t")[0]
        #id = eval(line.split("\t")[2])
        if name not in flag_list:
            out_list.append(line)
        else:
            continue
        flag_list.append(name)
    ouf = open(sys.argv[2], 'w')
    for i in range(len(out_list)):
        ouf.write(out_list[i] + "\n")
    inf.close()
    ouf.close()

try:
    main()
except IndexError:
    usage()

