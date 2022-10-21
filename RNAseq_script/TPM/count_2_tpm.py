#!/usr/bin/env python
'''
Copyright [yukaiquan 1962568272@qq.com]
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
    http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
input_file = `
Penicillin binding protein 3
Histidine decarboxylase
Nitric oxide synthase, inducible
'''

import sys
import numpy as np

# Read in the file
# gff_input: str = 'arahy.Tifrunner.gnm1.ann1.CCJH.gene_models_main.gff3'
gff_input: str = sys.argv[1]
# count_input: str = 'Arahy_all_tissue_htseq_count.csv'
count_input: str = sys.argv[2]
output: str = sys.argv[3]

gff_len = {}
with open(gff_input, 'r') as gff_file:
    for line in gff_file:
        if line.startswith('#'):
            continue
        else:
            line = line.strip().split('\t')
            if line[2] == 'gene':
                # 关键在提取gene长度,一般要用gene的真实长度比如transcript或者CDS但是因为是对gene进行定量所以就取了gene长度
                gene_id = line[8].split('Name=')[1].split(';')[0]
                gff_len[gene_id] = abs(int(line[4]) - int(line[3]))

count_matrix = np.loadtxt(count_input, delimiter=',', skiprows=1, dtype=str)
# count matrix
count_np = count_matrix[:, 1:].astype(np.float64)
count_index = count_matrix[:, 0]
count_length = np.array([gff_len[i] for i in count_index])
count_length = count_length.reshape((count_length.shape[0], 1))

# TPM matrix
RPK = count_np / (count_length / 1000)
# 深度标准化
PMSC_rpk = np.sum(RPK,axis=0)/1e6
PMSC_rpk = PMSC_rpk.reshape(PMSC_rpk.shape[0],1)
# 计算TPM
TPM_matrix = RPK/PMSC_rpk


name_matrix = np.concatenate((count_index.reshape((count_index.shape[0], 1)), TPM_matrix), axis=1)
name = np.loadtxt(count_input, delimiter=',', dtype=str)
name = name[0,:]
name_matrix = np.concatenate((name.reshape(1,name.shape[0]),name_matrix),axis=0)

np.savetxt(output, name_matrix, delimiter=',', fmt='%s')
