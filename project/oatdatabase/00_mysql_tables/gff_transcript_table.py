#!/use/bin/env python
import time
import sys

input_gff_list: list = ['ACD_sfs_sativa_sfs.gff3', 'CD_insularis_CN108634.gff3', 'A_longiglumis_CN58139.gff3', 'A_longiglumis_CN58138.gff3', 'A_strigosa_nc_sorted.gff3', 'A_atlantica_bmc_sorted.gff3',
                        'ACD_OT3098_v2_gene_annotations.gff3', 'ACD_sativa_sang.gff3', 'C_eriantha_bmc.gff3', 'CD_insularis_BYU209.gff3']

# output_oat: str = 'oat_transcript_table.txt'
output_oat: str = sys.argv[1]
# id_mapping: str = 'transcript_id_mapping.txt'
id_mapping: str = sys.argv[2]
'''
id,gene
1,A.satnudSFS1A01G000001
2,A.satnudSFS1A01G000002
3,A.satnudSFS1A01G000003
4,A.satnudSFS1A01G000004
5,A.satnudSFS1A01G000005
6,A.satnudSFS1A01G000006
7,A.satnudSFS1A01G000007
8,A.satnudSFS1A01G000008
'''
out_list: list = []
numpy_list: list = []
# 生成时间
time_info = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
# 生成基本table表
for input_gff in input_gff_list:
    with open(input_gff, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            else:
                line = line.strip().split('\t')
                if line[2] == 'mRNA':
                    # print(line[8].split('ID=')[1].split(';')[
                    #       0], line[0], line[3], line[4], line[8].split('Parent=')[1], sep='\t')
                    transcript_id = line[8].split('ID=')[1].split(';')[0]
                    chr = line[0]
                    start = str(line[3])
                    end = str(line[4])
                    strand = line[6]
                    mRNALength = str(abs(int(start) - int(end)))
                    if ';' in line[8].split('Parent=')[1]:
                        gene_id = line[8].split('Parent=')[1].split(';')[0]
                    else:
                        gene_id = line[8].split('Parent=')[1]
                    out_list.append(
                        [transcript_id, chr, start, end, strand, mRNALength, gene_id])
print('生成普通tab文件成功')
# 生成燕麦数据库中的基本表
id_mapping_dict: dict = {}
with open(id_mapping, 'r') as f:
    for line in f:
        line = line.strip().split(',')
        id_mapping_dict[line[1]] = line[0]
print('生成燕麦数据库专用mysql原始表')
with open(output_oat, 'w') as o:
    for i in out_list:
        if i[6] in id_mapping_dict:
            numpy_list.append([i[0], i[1], i[2], i[3], i[4], i[5],
                              '', time_info, time_info, id_mapping_dict[i[6]]])
            # o.write(i[0] + '\t' + i[1] + '\t' + i[2] + '\t' + i[3] + '\t' + i[4] + '\t' +
            #         '' + '\t' + time_info + '\t' + time_info + '\t' + id_mapping_dict[i[6]] + '\n')
        else:
            print('gene_id not in id_mapping_dict: ', i[6])
    for i in numpy_list:
        o.write(','.join(i) + '\n')

print('Done!')
