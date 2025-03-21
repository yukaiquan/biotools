#!/usr/bin/env python
# coding: utf-8
import sys
import os
import argparse
import subprocess


input_ipr_tsv = sys.argv[1]
# input_ipr_tsv = "./test_ipr.tsv"
input_go_basic_obo = sys.argv[2]
# input_go_basic_obo = "./go-basic.obo"
output_go_ann = sys.argv[3]
# output_go_ann = "test_ipr_go_ann.tsv"


go_ann_id_mapping = {}

with open(input_go_basic_obo, 'r') as f:
    go_id = ''
    go_name = ''
    namespace = ''
    def_desc = ''
    alt_ids = []
    for line in f:
        if line.startswith('id:'):
            go_id = line.split('id: ')[1].strip()
            # go_ann_id_mapping[go_id] = []
        elif line.startswith('name:'):
            go_name = line.split('name: ')[1].strip()
            # go_ann_id_mapping[go_id].append(go_name)
        elif line.startswith('namespace:'):
            namespace = line.split('namespace: ')[1].strip()
            # go_ann_id_mapping[go_id].append(namespace)
        elif line.startswith('def:'):
            def_desc = line.split('def: ')[1].strip().split('"')[1].strip()
            # go_ann_id_mapping[go_id].append(def_desc)
        elif line.startswith('alt_id:'):
            alt_id = line.split('alt_id: ')[1].strip()
            alt_ids.append(alt_id)
        elif line.startswith('[Term]'):
            if go_id != '':
                if go_id not in go_ann_id_mapping:
                    go_ann_id_mapping[go_id] = [go_name, namespace, def_desc]
                else:
                    print("go_idalready exists in go_ann_id_mapping: " + go_id)
                if len(alt_ids) > 0:
                    for alt_id in alt_ids:
                        go_ann_id_mapping[alt_id] = go_ann_id_mapping[go_id]
                go_id = ''
                go_name = ''
                namespace = ''
                def_desc = ''
                alt_ids = []
            else:
                continue
        else:
            continue

outpu_ann_dict = {}
with open(input_ipr_tsv, 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        else:
            lines = line.split('\t')
            gene_id = lines[0]
            pep_len = lines[2]
            align_start = lines[6]
            align_end = lines[7]
            e_value = lines[8]
            go_desc = lines[13]
            # print(go_desc)
            if "GO" in go_desc:
                if "|" in go_desc:
                    go_descs = go_desc.split("|")
                    for go_desc in go_descs:
                        if "GO" in go_desc:
                            go_id = go_desc.split("(")[0]
                            if go_id in go_ann_id_mapping:
                                go_ann = go_ann_id_mapping[go_id]
                                if gene_id not in outpu_ann_dict:
                                    outpu_ann_dict[gene_id] = [[go_id,go_ann,pep_len, align_start, align_end, e_value]]
                                else:
                                    outpu_ann_dict[gene_id].append([go_id,go_ann,pep_len, align_start, align_end, e_value])
                            else:
                                print("GO ID not found:", go_id)

                else:
                    go_id = go_desc.split("(")[0]
                    if go_id in go_ann_id_mapping:
                        go_ann = go_ann_id_mapping[go_id]
                        if gene_id not in outpu_ann_dict:
                            outpu_ann_dict[gene_id] = [[go_id,go_ann,pep_len, align_start, align_end, e_value]]
                        else:
                            outpu_ann_dict[gene_id].append([go_id,go_ann,pep_len, align_start, align_end, e_value])

with open(output_go_ann, 'w') as f:
    f.write("#GeneID\tGO_ID\tGO_Type\tGO_Annotation\tPeptide_Length\tAlign_Start\tAlign_End\tE_value\n")
    for gene_id in outpu_ann_dict.keys():
        # print(gene_id)
        go_anns = outpu_ann_dict[gene_id]
        # print(go_anns)
        for go_ann in go_anns:
            # print(go_ann)
            namespace = go_ann[1][1]
            ann = go_ann[1][2]
            f.write(gene_id + "\t" + go_ann[0] + "\t" + namespace + "\t" + ann + "\t" + go_ann[2] + "\t" + go_ann[3] + "\t" + go_ann[4] + "\t" + go_ann[5] + "\n")

    

    
