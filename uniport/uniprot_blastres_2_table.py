#!/usr/bin/env python
import sys
import gzip

# Read the UniProt BLAST output file
input_blast_out = sys.argv[2]
# Read the UniProt fasta
input_fasta = sys.argv[1]
output_file = sys.argv[3]

def read_fasta_2_ann(fasta_file)->dict:
    # 判断文件后缀是不是gz
    if fasta_file.endswith('.gz'):
        fasta_open = gzip.open(fasta_file, 'rt')
    else:
        fasta_open = open(fasta_file, 'r')
    
    fasta_dict = {}
    for line in fasta_open:
        line = line.strip()
        if line.startswith('>'):
            # 获取序列名
            seq_name = line.replace('>', '')
            seq_name = seq_name.split()[0]
            seq_ann = line.replace('>' + seq_name, '').strip()
            if seq_name not in fasta_dict:
                fasta_dict[seq_name] = seq_ann
            else:
                print('Warning: duplicate sequence name: ' + seq_name)
    fasta_open.close()
    return fasta_dict

fasta_dict = read_fasta_2_ann(input_fasta)
id_res = {}
with open(input_blast_out, 'r') as f:
    for line in f:
        line = line.strip()
        if line.startswith('#'):
            continue
        else:
            line_split = line.split('\t')
            query_id = line_split[0]
            subject_id = line_split[1]
            identical = float(line_split[2])
            query_len = int(line_split[3])
            evalue = float(line_split[10])
            if query_id not in id_res:
                content = fasta_dict[subject_id]
                id_res[query_id] = [content, identical,query_len, evalue,subject_id]
            else:
                print('Warning: duplicate query_id: ' + query_id)

with open(output_file, 'w') as f:
    f.write('query_id\tidentical\tquery_len\tevalue\tswissID\tannotation\n')
    for query_id in id_res:
        content = id_res[query_id][0]
        identical = id_res[query_id][1]
        query_len = id_res[query_id][2]
        evalue = id_res[query_id][3]
        subject_id = id_res[query_id][4]
        f.write(query_id + '\t' +str(identical) + '\t' + str(query_len) + '\t' + str(evalue) + '\t' +subject_id +'\t' + content + '\n')


