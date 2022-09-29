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
# input_file = r'D:\databasezip\genome\A\eriantha\chrUn.gff3'
output_file = sys.argv[2]
# output_file = r'D:\databasezip\genome\A\eriantha\test.gff3'

gff_input = np.loadtxt(input_file,
                       dtype='str', delimiter='\t')
# print(gff_input.size * gff_input.itemsize / 1e9, 'GB')
gff_df = pd.DataFrame(gff_input)
del gff_input
# sort gff by type and chr and position
gff_df[2] = gff_df[2].astype(cat_size_order)
gff_df[3] = gff_df[3].astype(int)
gff_df.sort_values(by=[0, 3, 4, 2], ascending=[
    True, True, True, True], inplace=True)
gff_df.reset_index(drop=True, inplace=True)
gff_df[3] = gff_df[3].astype(str)

rename_list = []
new_rename_list = []
exon_counts = {}
cds_counts = {}

# 提供哈希表因为有的基因存在跨基因的情况不能根据位置定义 exon 和 cds 的顺序
mRNA_cds_id_dict = {}
mRNA_exon_id_dict = {}

for line in gff_df.values:
    # line[8] = line[8].strip(';')
    if line[2] == 'gene':
        # 防止-RA -RB这种多个转录本乱入基因
        if '-' not in line[8]:
            # gene_id = line[8].split('=')[1]
            rename_list.append(line)
        else:
            print(line)
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
new_rename_list = np.array(new_rename_list)
np.savetxt(output_file,
           new_rename_list, fmt='%s', delimiter='\t')
