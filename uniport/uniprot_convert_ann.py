#! /use/bin/env python3
# -*- coding: utf-8 -*-
'''
chunk_size = 1024*1024*10
uniprot_dat file is a big file, so i want to read it in 10M size 
use python uniprot_convert_ann.py uniprot_sprot.dat uniprot_sprot.ann
'''
import sys


# input_file = "uniprot_plants_5000.dat"
# output_file = "uniprot_prot_ann.tsv"
input_file = sys.argv[1]
output_file = sys.argv[2]

""" 按块读取 """


def read_in_chunks(filePath: str, chunk_size: int = 1024*1024*1000):
    """
    Lazy function (generator) to read a file piece by piece.
    Default chunk size: 1024M
    You can set your own chunk size 
    """
    file_object = open(filePath)
    while True:
        chunk_data = file_object.read(chunk_size)
        if not chunk_data:
            break
        yield chunk_data

# uni_list = []


uniprot_id_list: list = []
for chunk in read_in_chunks(input_file):
    for line in chunk.split('\n'):
        if line.startswith("ID"):
            id = line.split('   ')[1]
        elif line.startswith("AC"):
            # print(line)
            count = 0
            if uniprot_id_list != []:
                with open(output_file, 'a+') as f:
                    # print(uniprot_id_list)
                    if len(uniprot_id_list) > 1:
                        for uni in uniprot_id_list:
                            f.write(uni + "\t" + id + "\t" + annotation + "\n")
                    else:
                        f.write(uniprot_id_list[0] + "\t" +
                                id + "\t" + annotation + "\n")
            if len(line.strip().split('   ')) > 1:
                uniprot_id_list = line.strip().split(
                    '   ')[1].replace(';', '').split()
                annotation: str = ""
            else:
                print('ID erro:', line)
        elif line.startswith("CC"):
            # print(line.split('   ')[1])
            # if 'P42652' in uniprot_id:
            #     print(line)
            if '--------' in line:
                continue
            elif 'Copyrighted by the UniProt Consortium' in line:
                count += 1
                if count > 1:
                    print(uniprot_id_list, count)
                continue
            elif 'Distributed under the Creative Commons Attribution' in line:
                continue
            else:
                line = line.replace('CC', '').replace('-!-', '').strip()
                # if 'P42652' in uniprot_id:
                #     print(line)
                annotation += line + " "
            # annotation += l.split('')[1].strip() + "\n"
        else:
            continue
