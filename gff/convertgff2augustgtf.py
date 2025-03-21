#!/usr/bin/rnv python
import sys
import re


# input_gff = sys.argv[1]
# # cds_fasta = sys.argv[2]
# output_gtf = sys.argv[2]

def convertgff2augustgtf(input_gff, output_gtf):
    transcript_dict = {}
    with open(input_gff, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith("#"):
                continue
            line_list = line.split("\t")
            if line_list[2] == "mRNA":
                transcript_id = line_list[8].split(";")[0].split("=")[1]
                if transcript_id not in transcript_dict:
                    gene_id = line_list[8].split(";")[1].split("=")[1]
                    transcript_dict[transcript_id] = gene_id
                else:
                    print("Error: duplicate transcript ID found: {}".format(transcript_id))
    # cds_code_dict = {}
    # cds_fasta_dict = {}
    # with open(cds_fasta, 'r') as f:
    #     for line in f:
    #         line = line.strip()
    #         if line.startswith(">"):
    #             cds_name = line.split()[0].split(">")[1]
    #             cds_fasta_dict[cds_name] = ""
    #         else:
    #             cds_fasta_dict[cds_name] += line

    # code_nucle_dict  = {}

    # for cds_name in cds_fasta_dict:
    #     start_code = ["ATG"]
    #     end_code = ["TAA", "TAG", "TGA"]
    #     cds = cds_fasta_dict[cds_name]
    #     for start_code_i in start_code:
    #         # 查找在cds中第一次出现start_code的位置使用正则表达式
    #         start_index = re.search(start_code_i, cds).start()
    #         # 查找在cds中最后一次出现end_code的位置使用正则表达式
    #         end_index = 0
    #         for end_code_i in end_code:
    #             end_index_s = re.search(end_code_i, cds[start_index:]).end()
    #             if end_index_s > end_index:
    #                 end_index = end_index_s
    #             else:
    #                 end_index = end_index
    #     code_nucle_dict[cds_name] = (start_index, end_index)

            
        
    out = []
    with open(input_gff, 'r') as f:
            start = 0
            end = 0
            strand = "+"
            cds_start = 0
            cds_end = 0
            intron_start = 0
            intron_end = 0
            for line in f:
                line = line.strip()
                # print(line)
                if line.startswith("#"):
                    out.append(line.split("\t"))
                else:
                    line_list = line.split("\t")
                    if line_list[2] == "gene":
                        line_list[8] = line_list[8].split("=")[1]
                        start = line_list[3]
                        end = line_list[4]
                        strand = line_list[6]
                        out.append(line_list)
                    elif line_list[2] == "mRNA":
                        transcript_id = line_list[8].split(";")[0].split("=")[1]
                        out.append([line_list[0], line_list[1], "transcript", line_list[3], line_list[4], line_list[5], line_list[6], line_list[7], transcript_id])
                        gene_id = transcript_dict[transcript_id]
                        line_list[4] = str(int(line_list[3]) + 2)
                        line_list[8] = 'transcript_id "{}"; gene_id "{}";'.format(transcript_id, gene_id)
                        if strand == "-":
                            line_list[2] = "stop_codon"
                        else:
                            line_list[2] = "start_codon"
                        out.append(line_list)
                    elif line_list[2] == "CDS" or line_list[2] == "exon":
                        transcript_id = line_list[8].split(";")[1].split("=")[1]
                        gene_id = transcript_dict[transcript_id]
                        line_list[8] = 'transcript_id "{}"; gene_id "{}";'.format(transcript_id, gene_id)
                        # out.write("\t".join(line_list) + "\n")
                        if line_list[2] == "exon":
                            start_exon = int(line_list[3])
                            end_exon = int(line_list[4])
                            if int(start_exon) > int(start) and int(end_exon) == int(end):
                                # 最后一个exon
                                out.append(line_list)
                                line_list[3] = str(int(line_list[4]) -2)
                                if strand == "+":
                                    line_list[2] = "stop_codon"
                                else:
                                    line_list[2] = "start_codon"
                                out.append(line_list)  
                            elif int(start_exon) == int(start) and int(end_exon) < int(end):
                                # 第一个exon 记录cds的终止位置为intron的起始位置
                                intron_start = end_exon + 1
                                out.append(line_list)
                                # line_list[2] = "intron"
                                # line_list[3] = str(intron_start)
                                # line_list[4] = str(start_exon - 1)
                                # out.write("\t".join(line_list) + "\n")
                            elif int(start_exon) == int(start) and int(end_exon) == int(end):
                                # 只有一个exon
                                out.append(line_list)
                                # line_list[3] = str(int(line_list[4]) -2)
                                # if strand == "+":
                                #     line_list[2] = "stop_codon"
                                # else:
                                #     line_list[2] = "start_codon"
                                # out.append(line_list)
                            else:
                                # 位于中间部分的exon
                                out.append(line_list)
                                line_list[2] = "intron"
                                line_list[3] = str(intron_start)
                                line_list[4] = str(start_exon - 1)
                                out.append(line_list)
                                intron_start = end_exon + 1
                        elif line_list[2] == "CDS":
                            start_cds = int(line_list[3])
                            end_cds = int(line_list[4])
                            # out.write("\t".join(line_list) + "\n")
                            if int(start_cds) > int(start) and int(end_cds) == int(end):
                                # 最后一个cds
                                type = "intron"
                                start_in = str(intron_start)
                                end_in = str(start_cds - 1)
                                # out.write(line_list[0] + "\t" + line_list[1] + "\t" + type + "\t" + start_in + "\t" + end_in + "\t" + strand + "\t.\t" + line_list[8] + "\n")
                                out.append([line_list[0], line_list[1], type, start_in, end_in,line_list[5], strand,line_list[7], line_list[8]])
                                out.append(line_list)
                            else:
                                out.append(line_list)
                    else:
                        print("Error: unexpected feature type found: {}".format(line_list[2]))
    

    # 去除数组中的重复
    unique_lines = unique_nested_arrays(out)
    # 对out进行排序 先按第0列和3列排序，第2列优先级gene > transcript > stop_code/start_code > CDS > exon
    # 对数据进行排序
    sorted_lines = sorted(unique_lines, key=sort_key)
    # 写入文件
    with open(output_gtf, 'w') as f:
        for line in sorted_lines:
            f.write("\t".join(line) + "\n")

# 定义排序的优先级
priority = {'gene': 1, 'transcript': 2, 'start_codon': 3, 'stop_codon': 3, 'CDS': 4, 'exon': 5}
# 使用自定义排序函数
def sort_key(line):
    # 获取染色体名称和起始位置
    chrom, start = line[0], int(line[3])
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


if __name__ == "__main__":
    try:
        input_gff = sys.argv[1]
        output_gtf = sys.argv[2]
        convertgff2augustgtf(input_gff, output_gtf)
    except Exception as e:
        print(e)
        print("Usage: python convertgff2augustgtf.py input_gff output_gtf")
        sys.exit(1)


