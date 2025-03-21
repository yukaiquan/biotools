#!/usr/bin/env python


import sys

input_pav_gene = sys.argv[1]

input_gfa = sys.argv[2]

output_gfa = sys.argv[3]

output_id_2_gene = sys.argv[4]

gene_pep_dict = {}

with open(input_pav_gene, 'r') as f:
    index = 1
    for line in f:
        if line.startswith('Gene'):
            continue
        else:
            gene_name = line.split('\t')[0]
            if gene_name not in gene_pep_dict:
                gene_pep_dict[gene_name] = 'PAN' + str(index).zfill(7)
                index += 1


with open(input_gfa, 'r') as f, open(output_gfa, 'w') as f_out:
    for line in f:
        lines = line.strip().split()
        length = len(lines)
        if length == 9:
            if lines[1] in gene_pep_dict:
                lines[1] = gene_pep_dict[lines[1]]
            name = lines[9].split(':')[2]
            if name in gene_pep_dict:
                taget_name = gene_pep_dict[name]
                lines[9] = lines[9].replace(name, taget_name)
            f_out.write('\t'.join(lines) + '\n')
        elif length == 11:
            if lines[1] in gene_pep_dict:
                lines[1] = gene_pep_dict[lines[1]]
            if lines[3] in gene_pep_dict:
                lines[3] = gene_pep_dict[lines[3]]
            f_out.write('\t'.join(lines) + '\n')
        else:
            if lines[0] == 'W':
                fasta_name = lines[6]
                # 以>或者<分割，都有
                
            else:
                f_out.write('\t'.join(lines) + '\n')


            


def replace_text_in_file(file_path,out_path, replacements):
    """
    根据字典内容替换文本文件中的内容。

    :param file_path: 文本文件的路径
    :param replacements: 一个字典，键是要被替换的字符串，值是替换后的字符串
    """
    # 读取原始文件内容
    with open(file_path, 'r', encoding='utf-8') as file:
        content = file.read()

    # 遍历字典中的每个键值对，进行替换
    for old_text, new_text in replacements.items():
        content = content.replace(old_text, new_text)

    # 将替换后的内容写回文件
    with open(out_path, 'w', encoding='utf-8') as file:
        file.write(content)


# replacements = {
#     'old_text1': 'new_text1',
#     'old_text2': 'new_text2',
#     'old_text3': 'new_text3'
# }

replace_text_in_file(input_gfa,output_gfa, gene_pep_dict)



