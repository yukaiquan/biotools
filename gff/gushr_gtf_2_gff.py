#!/usr/bin/env python

import sys
import re
# python /mnt/e/software/Ex-seq/biotools/gff/gushr_gtf_2_gff.py final_annotation.gff ACD_sativa_sfs_v2.gff3 ACD_sativa_sfs_v2_addutr.gff3

# gushr_gff = "D:\\databasezip\\genome\\ACD\\SFS\\09_addutr\\final_annotation.gff"
# old_gff = "D:\\databasezip\\genome\\ACD\\SFS\\09_addutr\\ACD_sativa_sfs_v2.gff3"
# out_gff = "D:\\databasezip\\genome\\ACD\\SFS\\09_addutr\\ACD_sativa_sfs_v2_addutr.gff3"
out_gff = sys.argv[3]
old_gff = sys.argv[2]
gushr_gff = sys.argv[1]


gff = {}
gene = {}
mrna = {}
cds = {}
cds_pos = {}
gene2mrna = {}
utr_5 = {}
utr_3 = {}

# 定义排序的优先级
priority = {'gene': 1, 'mRNA': 2, 'five_prime_UTR': 3, 'three_prime_UTR': 3, 'CDS': 4, 'exon': 5}
# 使用自定义排序函数
def sort_key(line):
    # 获取染色体名称和起始位置
    try:
        chrom, start = line[0], int(line[3])
    except:
        print(line)
    # 获取元件类型并根据优先级获取排序值
    feature = line[2]
    priority_value = priority.get(feature, 0)
    return (chrom, start, priority_value)

def unique_nested_arrays(nested_array):
    seen = set()
    unique = []

    for sub_array in nested_array:
        # 序列化子数组以便于比较
        tuple_sub_array = tuple(sub_array)
        if tuple_sub_array not in seen:
            seen.add(tuple_sub_array)
            unique.append(sub_array)

    return unique

def find_cds_num(name:str)->str:
    # text is CDS or exon
    match = re.search(r'(?<=[exon|CDS])\d',name)
    if match:
        return match.group(0)
    else:
        return ""

with open(gushr_gff,'r') as gff_file:
    for line in gff_file:
        if not line.startswith("#"):
            lines = line.strip().split("\t")
            if lines[2] == "gene":
                mrna_id = lines[8].split(";")[0].split("=")[1] 
                transcripts_num = lines[8].split(";")[1].split("=")[1]
                # 将最后的AVESA.00400a.r2.1Ag0000001.1的.1去掉 有多个.1 正则匹配
                pattern = r"\.[^.]+$"
                extracted_string = re.sub(pattern, '', mrna_id)
                lines[8] = f"ID={extracted_string}"
                lines[2] = "gene"
                if extracted_string in gene:
                    print(f"{extracted_string} has multiple gene")
                else:
                    gene[extracted_string] = lines
                if extracted_string not in gene2mrna:
                    gene2mrna[extracted_string] = [mrna_id]
                else:
                    gene2mrna[extracted_string].append(mrna_id)
                if int(transcripts_num) > 1:
                    print(f"{extracted_string} has {transcripts_num} transcripts")
            elif lines[2] == "prediction":
                mrna_id = lines[8].split(";")[0].split("=")[1].split("_")[0]
                pattern = r"\.[^.]+$"
                extracted_string = re.sub(pattern, '', mrna_id)
                lines[8] = f"ID={mrna_id};Parent={extracted_string}"
                lines[2] = "mRNA"
                if mrna_id in mrna:
                    print(f"{mrna_id} has multiple mrna")
                else:
                    mrna[mrna_id] = lines
            elif lines[2] == "five_prime_UTR" or lines[2] == "three_prime_UTR":
                mrna_name = lines[8].split("=")[1].split("_")[0]
                utr_type = ""
                strand = lines[6]
                if lines[2] == "five_prime_UTR":
                    utr_type = "utr5p"
                    lines[8] = f"ID={mrna_name}.{utr_type}.1;Parent={mrna_name}"
                    if mrna_name in utr_5:
                        # print(f"{mrna_name} has multiple utr5p")
                        start_utr = utr_5[mrna_name][3]
                        end_utr = utr_5[mrna_name][4]
                        if strand == "-":
                            if int(lines[3]) < int(start_utr):
                                # utr_5[mrna_name] = (lines[3], end_utr)
                                lines[4] = end_utr
                            else:
                                print(f"{mrna_name} utr5p is error!")
                        else:
                            if int(lines[4]) > int(end_utr):
                                # utr_5[mrna_name] = (start_utr, lines[4])
                                lines[3] = start_utr
                            else:
                                print(f"{mrna_name} utr5p is error!")
                        utr_5[mrna_name] = lines
                    else:
                        # utr_5[mrna_name] = (lines[3], lines[4])
                        utr_5[mrna_name] = lines
                elif lines[2] == "three_prime_UTR":
                    utr_type = "utr3p"
                    lines[8] = f"ID={mrna_name}.{utr_type}.1;Parent={mrna_name}"
                    if mrna_name in utr_3:
                        # print(f"{mrna_name} has multiple utr3p")
                        start_utr = utr_3[mrna_name][3]
                        end_utr = utr_3[mrna_name][4]
                        if strand == "-":
                            if int(lines[3]) < int(start_utr):
                                # utr_3[mrna_name] = (lines[3], end_utr)
                                lines[4] = end_utr
                            else:
                                print(f"{mrna_name} utr3p is error!")
                                # utr_3[mrna_name] = (start_utr, lines[4])
                        else:
                            if int(lines[4]) > int(end_utr):
                                # utr_3[mrna_name] = (start_utr, lines[4])
                                lines[3] = start_utr
                            else:
                                print(f"{mrna_name} utr3p is error!")
                        utr_3[mrna_name] = lines
                    else:
                        # utr_3[mrna_name] = (lines[3], lines[4])
                        utr_3[mrna_name] = lines
                # lines[8] = f"ID={mrna_name}.{utr_type}.1;Parent={mrna_name}"
                # if mrna_name in utr:
                #     utr[mrna_name].append(lines)
                # else:
                #     utr[mrna_name] = [lines]
            elif lines[2] == "CDS":
                mrna_name = lines[8].split(";")[1].split("=")[1].split("_")[0]
                cds_id = lines[8].split(";")[0].split("=")[1]
                cds_id = cds_id.replace("CDS","CDS.")

                lines[8] = f"ID={cds_id};Parent={mrna_name}"
                if mrna_name in gff:
                    cds[mrna_name].append(lines)
                else:
                    cds[mrna_name] = [lines]
                chr = lines[0]
                start = lines[3]
                end = lines[4]
                name = chr + ":" + start + "-" + end
                if name in cds_pos:
                    print(f"{name} has multiple cds")
                else:
                    cds_pos[name] = lines

res = []

with open(old_gff,'r') as old_gff_file:
    mrna_start =  ""
    mrna_end = ""
    old_mrn_start = ""
    old_mrn_end = ""
    for line in old_gff_file:
        if not line.startswith("#"):
            lines = line.strip().split("\t")
            if lines[2] == "gene":
                gene_id = lines[8].split("=")[1]
                if gene_id in gene:
                    gene_line = gene[gene_id]
                    # out_gff_file.write("\t".join(gene_line)+"\n")
                    res.append(gene_line)
                else:
                    print(f"{gene_id} not in gushr")
                    res.append(lines)
            elif lines[2] == "mRNA":
                mrna_id = lines[8].split(";")[0].split("=")[1]
                old_mrn_start = lines[3]
                old_mrn_end = lines[4]
                if mrna_id in mrna:
                    mrna_line = mrna[mrna_id]
                    # out_gff_file.write("\t".join(mrna_line)+"\n")
                    res.append(mrna_line)
                    mrna_start = mrna_line[3]
                    mrna_end = mrna_line[4]
                else:
                    print(f"{mrna_id} not in gushr")
                    res.append(lines)
                    mrna_start = lines[3]
                    mrna_end = lines[4]
            elif lines[2] == "CDS":
                res.append(lines)
            elif lines[2] == "exon":
                cds_id = lines[8].split(";")[1].split("=")[1]
                if mrna_start == "":
                    print(f"{cds_id} not in mrna, and is error!")
                else:
                    if int(old_mrn_start) == int(lines[3]):
                        lines[3] = mrna_start
                    elif int(old_mrn_end) == int(lines[4]):
                        lines[4] = mrna_end
                res.append(lines)
            else:
                print(f"{lines[2]} is unknower type!")

for i in res:
    if i[2] == "CDS":
        cds_id = i[8].split(";")[0].split("=")[1]
        cds_num = find_cds_num(cds_id)
        new_cds_id = cds_id.replace("CDS" + cds_num + ".", "")
        new_cds_id = new_cds_id + ".CDS." + cds_num
        i[8] = i[8].replace(cds_id, new_cds_id)
    elif i[2] == "exon":
        cds_id = i[8].split(";")[0].split("=")[1]
        cds_num = find_cds_num(cds_id)
        new_cds_id = cds_id.replace("exon" + cds_num + ".", "")
        new_cds_id = new_cds_id + ".exon." + cds_num
        i[8] = i[8].replace(cds_id, new_cds_id)


for u in utr_5:
    utr_line = utr_5[u]
    res.append(utr_line)

for u in utr_3:
    utr_line = utr_3[u]
    res.append(utr_line)
                
sorted_lines = sorted(res, key=sort_key)

with open(out_gff, 'w') as f:
    for line in sorted_lines:
        f.write("\t".join(line) + "\n")
                










            