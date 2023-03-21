#!/usr/bin/env python
# 读取压缩文件
import gzip
from tqdm import tqdm
import sys

# python ../03_vcf/snpEff_ann_2_table.py ../03_vcf/ogle.filter.ann.vcf.gz ogle.filter.ann2.txt oat1ogle

# input_file = "ogle.filter.ann.vcf.gz"
input_file = sys.argv[1]
# output_file = "ogle.filter.ann.vcf.gz.table"
output_file = sys.argv[2]
prefix = sys.argv[3]


# 压缩文件太大，需要分批读取
# 分批读取压缩文件
def read_gzip_file(input_file, batch_size=100000):
    with gzip.open(input_file, 'rb') as f:
        while True:
            lines = f.readlines(batch_size)
            if not lines:
                break
            for line in lines:
                yield line


input_list = []

# 读取压缩文件
for line in tqdm(read_gzip_file(input_file), desc='Reading file'):
    line = line.decode('utf-8')
    if line.startswith('#'):
        continue
    cols = line.split('\t')
    # 第六列为PASS保留
    if cols[6] == 'PASS':
        chr_name = cols[0]
        pos = cols[1]
        marker_name = prefix + '_' + \
            chr_name.replace('chr', '') + '_' + str(pos)
        ref = cols[3]
        alt = cols[4]
        # 判断变异类型
        # if len(ref) == len(alt):
        #     if len(ref) == 1:
        #         snp_type = 'snp'
        #     else:
        #         snp_type = 'mnp'
        # elif len(ref) > len(alt):
        #     snp_type = 'del'
        # elif ',' in alt:
        #     snp_type = 'snp'
        # else:
        #     snp_type = 'ins'
        if ',' in alt:
            # 杂合突变
            snp_type = 'het'
        elif len(ref) == len(alt):
            if len(ref) == 1:
                snp_type = 'snp'
            else:
                snp_type = 'mnp'
        elif len(ref) > len(alt):
            snp_type = 'del'
        else:
            snp_type = 'ins'
        qual = cols[5]
        # 第二列为snpEff的注释信息
        # AC=2;AF=1;AN=2;DP=2;ExcessHet=3.0103;FS=0;MLEAC=1;MLEAF=0.5;MQ=40;QD=31.63;SOR=2.303;ANN=A|downstream_gene_variant|MODIFIER|A.satnudSFS1A01G007012|A.satnudSFS1A01G007012|transcript|A.satnudSFS1A01G007012.1|protein_coding||c.*2050G>T|||||2050|,A|intergenic_region|MODIFIER|CHR_START-A.satnudSFS1A01G007012|CHR_START-A.satnudSFS1A01G007012|intergenic_region|CHR_START-A.satnudSFS1A01G007012|||n.249948C>A||||||
        # 逐个注释信息操作
        ann = cols[7].split(';')
        # 提取QD、SOR、FS、MQ、BaseQRankSum、ReadPosRankSum 没有的用NULL填充
        dp = 'NULL'
        qd = 'NULL'
        sor = 'NULL'
        fs = 'NULL'
        mq = 'NULL'
        baseqranksum = 'NULL'
        readposranksum = 'NULL'
        # 位置不确定，需要遍历
        for i in ann:
            # print(i)
            # if i.startswith('QD'): 识别不到QD
            if i.startswith('QD='):
                # print(i)
                qd = i.split('=')[1]
                # print(qd)
            elif i.startswith('DP='):
                dp = i.split('=')[1]
            elif i.startswith('SOR='):
                sor = i.split('=')[1]
            elif i.startswith('FS='):
                fs = i.split('=')[1]
            elif i.startswith('MQ='):
                mq = i.split('=')[1]
            elif i.startswith('BaseQRankSum='):
                baseqranksum = i.split('=')[1]
            elif i.startswith('ReadPosRankSum='):
                readposranksum = i.split('=')[1]
            else:
                continue
        ann_eff = ann[-1].split('|')
        # 取出注释信息中的第二三个字段
        ann_pos = ann_eff[1]
        ann_effect = ann_eff[2]
        # 写入列表
        input_list.append([marker_name, chr_name, pos, ref, alt, snp_type, qual,
                          dp, qd, sor, fs, mq, baseqranksum, readposranksum, ann_pos, ann_effect])

# 写入文件
with open(output_file, 'w') as f:
    # 写入标题
    f.write('marker_name\tchr\tpos\tref\talt\tsnp_type\tqual\tdp\tqd\tsor\tfs\tmq\tbaseqranksum\treadposranksum\tsnp_gene_pos\tsnp_effect\n')
    # 逐行写入
    for i in tqdm(input_list, desc='Writing file'):
        f.write('\t'.join(i) + '\n')
