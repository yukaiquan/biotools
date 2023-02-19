#!/usr/bin/env python3
import sys
import os
import time

# 该代码的功能是将orthofind输出的结果中的一对多关系提取出来并且将同一物种的一对多关系合并在一起

# input_dir = "./Orthologues_ACD.SFS.pep/"
input_dir = sys.argv[1]
# input_dir = '/home/zhao/orthofinder/Orthologues_ACD.SFS.pep/'
output_file = sys.argv[2]
# 需要排除的文件列表
exclude_files = sys.argv[3].split(",")
# 判断exclude_files是否为空
if exclude_files[0] == "":
    exclude_files = []


def get_orthologs(ortholog_file: str) -> dict:
    orthologs = {}
    with open(ortholog_file, "r") as f:
        for line in f:
            if line.startswith("Orthogroup"):
                continue
            line = line.strip().split("\t")
            if ',' in line[1]:
                for i in line[1].split(","):
                    if i not in orthologs:
                        orthologs[i.strip()] = line[2].replace(" ", "")
                    else:
                        print("Error: duplicate gene ID", i)
            else:
                if line[1] not in orthologs:
                    orthologs[line[1].strip()] = line[2].strip().replace(
                        " ", "")
                else:
                    print("Error: duplicate gene ID", line[1])
    return orthologs


# 从orthofind输出结果中提取一对多的关系
#orthologs: dict = get_orthologs(input_file)
# 依次读取orthofind输出结果中的每个文件，提取一对多的关系
# 列出文件夹中的所有文件
files = os.listdir(input_dir)
# 筛选出tsv文件
tsv_files = [i for i in files if i.endswith('.tsv')]
# 排除exclude_files中的文件
tsv_files = [input_dir + '/' + i for i in tsv_files if i not in exclude_files]
# 依次读取tsv文件
orthologs_all = {}
for tsv_file in tsv_files:
    orthologs: dict = get_orthologs(tsv_file)
    # 将每个tsv文件中的一对多关系添加到orthologs_all中
    # 如果orthologs_all中已经存在该键，则更新值
    for key, value in orthologs.items():
        if key not in orthologs_all:
            orthologs_all[key] = value
        else:
            orthologs_all[key] = orthologs_all[key] + ',' + value
# 获取当前时间
now = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())

# 按照键对orthologs_all进行排序
orthologs_all = dict(sorted(orthologs_all.items(), key=lambda x: x[0]))
with open(output_file, "w") as f:
    for key, value in orthologs_all.items():
        f.write(key + "\t" + value + "\t" + now + "\t" + now + "\n")
