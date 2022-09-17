'''
Author: yukaiquan
Date: 2022-02-13 18:50:06
Email: 1962568272@qq.com
LastEditor: yukaiquan
Description: Do not edit
'''
# -*- coding : utf-8 -*-

import argparse
import os
import re
import time


def dirList(dir):
    t = open(dir, "r")
    filelist = t.read().splitlines()
    t.close()
    allfile = []
    for filename in filelist:
        with open(filename, 'r') as file:
            lines = file.readlines()  # 读取每一行
        a = ''  # 空字符（中间不加空格）
        e = ''
        for line in lines:
            if line[0] == '>':
                e += line
            else:
                a += line.strip()  # strip()是去掉每行末尾的换行符\n 1
        c = a.split()  # 将a分割成每个字符串 2
        b = ''.join(c)  # 将c的每个字符不以任何符号直接连接 3
        with open(filename, 'w') as file:
            file.write(e + b)
        #以上部分将碱基序列合并为一行#
        f = open(filename)
        matchPattern = re.compile(r'>')
        while 1:
            line = f.readline()
            if not line:
                #print("Read file End or Error")
                break
            elif matchPattern.search(line):
                pass
            else:
                allfile.append(line)
        f.close()
    # allfile.append(f.read())
    return allfile


def Nfa(nnumbers, allfile):
    A = nnumbers * 'N'
    pathjar = A.join('%s' % a for a in allfile)
    # '%s' %a for a in
    list = []
    for line in pathjar:
        line = line.strip('\n')
        list.append(line)
    seqN = ''.join(list)
    result1 = [seqN[i:i + 60] for i in range(0, len(seqN), 60)]
    result = '\n'.join(result1)
    return result


def rename(filename, result):
    with open(filename + '.fasta', 'w') as file:
        file.write('>' + filename + '\n' + result)


def main():
    parser = argparse.ArgumentParser(
        description='此脚本用于链接多条序列并以多个N分割 最终再以默认一行60个显示的fasta文件')
    parser.add_argument(
        '-l', '--list', help='输入待合并的文件名列表注意：文件名需要全称', type=str, required=True)
    parser.add_argument('-n', '--N_number', help='输入两个序列间N的长度,默认以100N为分割',
                        type=int, required=False, default=100)
    parser.add_argument('-f', '--fastaname', help='输出fasta文件的名字 如果不提供将默认chrUn',
                        type=str, required=False, default='chrUn')
    args = parser.parse_args()
    time_start = time.time()
    print('开始按顺序批量合并文件')
    allfile = dirList(dir=args.list)
    print('开始为每个文件末尾添加N并按每行60字符串换行')
    result = Nfa(nnumbers=args.N_number, allfile=allfile)
    print('为新序列添加名字')
    rename(filename=args.fastaname, result=result)
    time_end = time.time()
    print('运行结束！========time cost', time_end - time_start, 's')


if __name__ == "__main__":
    main()
