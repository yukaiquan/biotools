#!/usr/bin/env bash
import sys
import os

fasta_kmer = sys.argv[1]
input_gz_vcf = sys.argv[2]

kmer_dict = {}
with open(fasta_kmer, 'r') as f:
    kmer_id = None
    kmer_seq = None
    for line in f:
        if line.startswith('>'):
            kmer_id = line.strip().lstrip('>')
        else:
            kmer_seq = line.strip()
            if kmer_id in kmer_dict:
                print(f"Warning: Duplicate kmer ID found: {kmer_id}", file=sys.stderr)
                sys.exit(1)
            kmer_dict[kmer_id] = kmer_seq

import gzip
res_list = []
res_dict = {}
with gzip.open(input_gz_vcf, 'rt') as f:
    header = []
    for line in f:
        if line.startswith('#'):
            print(line.strip())
            # 查找#CHROM开头的行
            if line.startswith('#CHROM'):
                header = line.strip().split('\t')
            continue
        fields = line.strip().split('\t')
        chr_name = fields[0]
        pos = fields[1]
        ref = fields[3]
        alt = fields[4]
        samples_length = len(fields) - 9
        # chr1A:1686299:G
        kmer_ref_id = f"{chr_name}:{pos}:{ref}"
        kmer_alt_id = f"{chr_name}:{pos}:{alt}"
        if kmer_ref_id not in kmer_dict:
            print(f"Error: Reference kmer ID {kmer_ref_id} not found in FASTA file", file=sys.stderr)
            # sys.exit(1)
            continue
        elif kmer_alt_id not in kmer_dict:
            print(f"Error: Alternate kmer ID {kmer_alt_id} not found in FASTA file", file=sys.stderr)
            # sys.exit(1)
            continue
        else:
            kmer_ref_seq = kmer_dict[kmer_ref_id]
            kmer_alt_seq = kmer_dict[kmer_alt_id]
            for i in range(9, 9 + samples_length):
                genotype = fields[i].split(':')[0]
                dp_ref = fields[i].split(':')[1] if len(fields[i].split(':')) > 1 else '.'
                dp_alt = fields[i].split(':')[2] if len(fields[i].split(':')) > 2 else '.'
                # 如果dp是.，则替换为0
                if dp_ref == '.':
                    dp_ref = '0'
                if dp_alt == '.':
                    dp_alt = '0'
                sample_name = header[i]
                if genotype == '0/0':
                    # res_list.append((sample_name, chr_name,pos, dp_ref,"0"))
                    if sample_name not in res_dict:
                        res_dict[sample_name] = [(chr_name,pos, dp_ref,"0")]
                    else:
                        res_dict[sample_name].append((chr_name,pos, dp_ref,"0"))
                elif genotype == '0/1':
                    # res_list.append((sample_name, chr_name,pos, dp_ref,dp_alt))
                    if sample_name not in res_dict:
                        res_dict[sample_name] = [(chr_name,pos, dp_ref,dp_alt)]
                    else:
                        res_dict[sample_name].append((chr_name,pos, dp_ref,dp_alt))
                elif genotype == '1/1':
                    # res_list.append((sample_name, chr_name,pos, "0",dp_alt))
                    if sample_name not in res_dict:
                        res_dict[sample_name] = [(chr_name,pos, "0",dp_alt)]
                    else:
                        res_dict[sample_name].append((chr_name,pos, "0",dp_alt))
                else:
                    # print(f"Warning: Unrecognized genotype {genotype} for sample {sample_name} at {chr_name}:{pos}", file=sys.stderr)
                    continue

# 输出结果
for sample in res_dict.keys():
    # 新建文件夹
    os.makedirs(sample, exist_ok=True)
    # 压缩写入
    samples_file = os.path.join(sample, f"{sample}_gbs_count.tsv.gz")
    with gzip.open(samples_file, 'wt') as out_f:
        for record in res_dict[sample]:
            out_f.write(f"{record[0]}\t{record[1]}\t{record[2]}\t{record[3]}\n")