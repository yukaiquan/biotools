#!/usr/bin/env python


import sys

input_file_samples = sys.argv[1]

out_put_fasta = sys.argv[2]

out_id_2_gene = sys.argv[3]

# AVENA.00001a.r1.PAN0000001

gene_pep_dict = {}

with open(input_file_samples, 'r') as f:
    for line in f:
        if line.startswith('>'):
            gene_id = line[1:]
            gene_pep_dict[gene_id] = ''
        else:
            gene_pep_dict[gene_id] += line

out_id_open = open(out_id_2_gene, 'w')

with open(out_put_fasta, 'w') as f:
    index = 1
    for gene_id in gene_pep_dict:
        name = 'AVENA.00001a.r1.PAN' + str(index).zfill(7)
        out_id_open.write(name+'\t'+gene_id+'\n')
        f.write('>'+name+'\n')
        f.write(gene_pep_dict[gene_id])
        index += 1

out_id_open.close()