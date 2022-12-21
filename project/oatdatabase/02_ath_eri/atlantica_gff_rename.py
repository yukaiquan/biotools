#! /usr/bin/env python3
'''
20220920 yukaiquan 1962568272@qq.com
chr1A   CoGe    gene    548208  548907  .   -   .   ID=AA006984
chr1A   CoGe    mRNA    548208  548907  .   -   .   ID=AA006984-RA;Parent=AA006984
chr1A   CoGe    CDS 548373  548486  .   -   1   ID=AA006984-RA:CDS;Parent=AA006984-RA
chr1A   CoGe    exon    548373  548486  .   -   .   ID=AA006984-RA:exon;Parent=AA006984-RA
chr1A   CoGe    CDS 548615  548736  .   -   0   ID=AA006984-RA:CDS;Parent=AA006984-RA
chr1A   CoGe    exon    548615  548736  .   -   .   ID=AA006984-RA:exon;Parent=AA006984-RA
chr1A   CoGe    CDS 548833  548907  .   -   0   ID=AA006984-RA:CDS;Parent=AA006984-RA
chr1A   CoGe    exon    548833  548907  .   -   .   ID=AA006984-RA:exon;Parent=AA006984-RA
chr1A   CoGe    CDS 578850  578951  .   -   0   ID=AA006985-RA:CDS;Parent=AA006985-RA
chr1A   CoGe    exon    578850  578951  .   -   .   ID=AA006985-RA:exon;Parent=AA006985-RA

I will try to rename CDS/exon numbers 
eg:
chr1A   CoGe    gene    548208  548907  .   -   .   ID=AA006984
chr1A   CoGe    mRNA    548208  548907  .   -   .   ID=AA006984-RA;Parent=AA006984
chr1A   CoGe    exon    548208  548253  .   -   .   ID=AA006984-RA:exon.1;Parent=AA006984-RA
chr1A   CoGe    CDS 548208  548253  .   -   1   ID=AA006984-RA:CDS.1;Parent=AA006984-RA
chr1A   CoGe    exon    548373  548486  .   -   .   ID=AA006984-RA:exon.2;Parent=AA006984-RA
chr1A   CoGe    CDS 548373  548486  .   -   1   ID=AA006984-RA:CDS.2;Parent=AA006984-RA
chr1A   CoGe    exon    548615  548736  .   -   .   ID=AA006984-RA:exon.3;Parent=AA006984-RA
chr1A   CoGe    CDS 548615  548736  .   -   0   ID=AA006984-RA:CDS.3;Parent=AA006984-RA
chr1A   CoGe    exon    548833  548907  .   -   .   ID=AA006984-RA:exon.4;Parent=AA006984-RA
chr1A   CoGe    CDS 548833  548907  .   -   0   ID=AA006984-RA:CDS.4;Parent=AA006984-RA


'''


import numpy as np
import pandas as pd
from pandas.api.types import CategoricalDtype
import sys

# define the order of the categories
cat_size_order = CategoricalDtype(
    ['gene', 'mRNA', 'exon', 'CDS'],
    ordered=True
)
input_file = sys.argv[1]
# input_file = r'D:\databasezip\genome\A\eriantha\genome20221214\new.gff3.gz'
output_file = sys.argv[2]
# output_file = r'D:\databasezip\genome\A\eriantha\genome20221214\test.gff3'

gff_input = np.loadtxt(input_file,
                       dtype='str', delimiter='\t')
# print(gff_input.size * gff_input.itemsize / 1e9, 'GB')
gff_df = pd.DataFrame(gff_input)
del gff_input
# sort gff by type and chr and position
gff_df[2] = gff_df[2].astype(cat_size_order)
gff_df[3] = gff_df[3].astype(int)
gff_df.sort_values(by=[0, 3, 2], ascending=[
    True, True, True], inplace=True)
gff_df.reset_index(drop=True, inplace=True)
gff_df[3] = gff_df[3].astype(str)

rename_list = []
new_rename_list = []
exon_counts = {}
cds_counts = {}

# 提供哈希表因为有的基因存在跨基因的情况不能根据位置定义 exon 和 cds 的顺序
mRNA_cds_id_dict = {}
mRNA_exon_id_dict = {}
# 准备一个基因dict 存储数据存在多个基因情况我们需要最长基因
gene_dict = {}
for line in gff_df.values:
    if line[2] == 'gene':
        if '-' in line[8].split('=')[1].replace(';', ''):  # 防止-RA -RB这种转录本编号乱入基因
            line[8] = line[8].split('-')[0] + ';'
            gene_id = line[8].split('=')[1].replace(';', '').split('-')[0]
        else:
            gene_id = line[8].split('=')[1].replace(';', '')
        if gene_id not in gene_dict:
            gene_dict[gene_id] = line
        else:
            if abs(int(line[4]) - int(line[3])) > abs(int(gene_dict[gene_id][4]) - int(gene_dict[gene_id][3])):
                gene_dict[gene_id] = line
            else:
                print('gene_id', gene_id, 'is not the longest gene')
    else:
        continue

for line in gff_df.values:
    # line[8] = line[8].strip(';')
    if line[2] == 'gene':
        # 防止-RA -RB这种多个转录本乱入基因
        if line[8].split('=')[1].replace(';', '') in gene_dict:
            rename_list.append(
                gene_dict[line[8].split('=')[1].replace(';', '')])
        else:
            print(line[8], 'is not gene')
    elif line[2] == 'mRNA':
        mRNA_id = line[8].split(';')[0].split('=')[1]
        mRNA_cds_id_dict[mRNA_id] = 1
        mRNA_exon_id_dict[mRNA_id] = 1
        rename_list.append(line)
    elif line[2] == 'exon':
        # exon_id = line[8].split(';')[0].split('=')[1]
        parent_id = line[8].split(';')[1].split('=')[1]

        exon_counts[parent_id] = exon_counts.get(parent_id, 0) + 1
        rename_list.append(line)
    elif line[2] == 'CDS':
        # CDS_id = line[8].split(';')[0].split('=')[1]
        parent_id = line[8].split(';')[1].split('=')[1]
        cds_counts[parent_id] = cds_counts.get(parent_id, 0) + 1
        rename_list.append(line)
    else:
        print('Error: unknown type')
        print(line)
        rename_list.append(line)
del gff_df


for line in rename_list:
    if line[2] == 'gene':
        # 防止-RA -RB这种多个转录本乱入基因
        if '-' not in line[8]:
            # if line[8].split('=')[1] == 'AE040502':
            #     print(line)
            gene_id = line[8].split('=')[1]
            # gene_id = line[8].split('=')[1]
            new_rename_list.append(line)
        else:
            print(line)
    elif line[2] == 'mRNA':
        # if line[8].split(';')[1].split('=')[1] == 'AE040502':
        #     print(line)
        mRNA_id = line[8].split(';')[0].split('=')[1]
        new_rename_list.append(line)
    elif line[2] == 'exon':
        parent_id = line[8].split(';')[1].split('=')[1]
        # if parent_id == 'AE040502-RA':
        #     print(mRNA_exon_id_dict[parent_id], exon_counts[parent_id])
        #     print(line)
        if exon_counts[parent_id] == 1:
            new_rename_list.append(line)
        else:
            line[8] = line[8].replace(
                'exon', 'exon' + '.' + str(mRNA_exon_id_dict[parent_id]))
            mRNA_exon_id_dict[parent_id] += 1
            new_rename_list.append(line)
    elif line[2] == 'CDS':
        parent_id = line[8].split(';')[1].split('=')[1]
        # if parent_id == 'AE040502-RA':
        #     print(mRNA_cds_id_dict[parent_id], cds_counts[parent_id])
        #     print(line)
        if cds_counts[parent_id] == 1:
            new_rename_list.append(line)
        else:
            line[8] = line[8].replace(
                'CDS', 'CDS' + '.' + str(mRNA_cds_id_dict[parent_id]))
            mRNA_cds_id_dict[parent_id] += 1
            new_rename_list.append(line)
    else:
        new_rename_list.append(line)

# 去除转录本中的-RA -RB
for line in new_rename_list:
    if line[2] == "mRNA":
        if '-' in line[8].split(";")[1]:
            line[8] = line[8].split(";")[0] + ";" + \
                line[8].split(";")[1].split('-')[0] + ";"
        else:
            line[8] = line[8]
# 再次进行排序
gff_df = pd.DataFrame(new_rename_list)
# sort gff by type and chr and position
gff_df[2] = gff_df[2].astype(cat_size_order)
gff_df[3] = gff_df[3].astype(int)
gff_df.sort_values(by=[0, 3, 2], ascending=[
    True, True, True], inplace=True)
gff_df.reset_index(drop=True, inplace=True)
gff_df[3] = gff_df[3].astype(str)

# 写入文件
gff_df.to_csv(output_file, sep='\t', header=False, index=False)
