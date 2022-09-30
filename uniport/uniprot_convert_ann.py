#! /use/bin/env python3
# -*- coding: utf-8 -*-
import sys
'''
chunk_size = 1024*1024*10
uniprot_dat file is a big file, so i want to read it in 10M size 
use python uniprot_convert_ann.py uniprot_sprot.dat uniprot_sprot.ann
'''


# input_file = "uniprot_plants_5000.dat"
# output_file = "uniprot_prot_ann.tsv"
input_file = sys.argv[1]
output_file = sys.argv[2]

""" 按块读取 """


def read_in_chunks(filePath, chunk_size=1024*1024*10):
    """
    Lazy function (generator) to read a file piece by piece.
    Default chunk size: 1M
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
        if line.startswith("AC"):
            # print(line)
            count = 0
            if uniprot_id_list != []:
                with open(output_file, 'a+') as f:
                    # print(uniprot_id_list)
                    if len(uniprot_id_list) > 1:
                        for uni in uniprot_id_list:
                            f.write(uni + "\t" + annotation + "\n")
                    else:
                        f.write(uniprot_id_list[0] + "\t" + annotation + "\n")
            uniprot_id = line.strip().split('   ')[1].split(';')
            if '' in uniprot_id:
                uniprot_id.remove('')
            uniprot_id_list = list(map(lambda x: x.strip(), uniprot_id))
            # uniprot_id_str = ','.join(uniprot_id)
            annotation: str = ""
        elif line.startswith("CC"):
            # print(line.split('   ')[1])
            # if 'P42652' in uniprot_id:
            #     print(line)
            if '--------' in line:
                continue
            elif 'Copyrighted by the UniProt Consortium' in line:
                count += 1
                if count > 1:
                    print(uniprot_id, count)
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
        # print(annotation)


# with open(input_file, 'r') as f:
#     lines = f.readlines()
#     # print(lines)
#     # lines = [line.decode('utf-8') for line in lines]
#     # chunks = read_file_by_chunk(input_file)
#     uniprot_id_list: list = []
#     for line in lines:
#         if line.startswith("AC"):
#             # print(line)
#             count = 0
#             if uniprot_id_list != []:
#                 with open(output_file, 'w') as f:
#                     if len(uniprot_id) > 1:
#                         for uni in uniprot_id:
#                             f.write(uni + "\t" + annotation + "\n")
#                     else:
#                         f.write(uniprot_id[0] + "\t" + annotation + "\n")
#             uniprot_id = line.strip().split('   ')[1].split(';')
#             if '' in uniprot_id:
#                 uniprot_id.remove('')
#             uniprot_id_list = list(map(lambda x: x.strip(), uniprot_id))
#             # uniprot_id_str = ','.join(uniprot_id)
#             annotation: str = ""
#         elif line.startswith("CC"):
#             # print(line.split('   ')[1])
#             # if 'P42652' in uniprot_id:
#             #     print(line)
#             if '--------' in line:
#                 continue
#             elif 'Copyrighted by the UniProt Consortium' in line:
#                 count += 1
#                 if count > 1:
#                     print(uniprot_id, count)
#                 continue
#             elif 'Distributed under the Creative Commons Attribution' in line:
#                 continue
#             else:
#                 line = line.replace('CC', '').replace('-!-', '').strip()
#                 # if 'P42652' in uniprot_id:
#                 #     print(line)
#                 annotation += line + " "
#             # annotation += l.split('')[1].strip() + "\n"
#         else:
#             continue
#         # print(annotation)


# with open(output_file, 'w') as f:
#     for uni in uni_list:
#         if len(uni[0]) == 1:
#             f.write(uni[0][0] + "\t" + uni[1] + "\n")
#         else:
#             for u in uni[0].split(','):
#                 f.write(u + "\t" + uni[1] + "\n")
