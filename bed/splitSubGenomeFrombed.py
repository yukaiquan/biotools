#!/usr/bin/env python

import sys
import re

input_bed = sys.argv[1]
input_cds = sys.argv[2]
input_pep = sys.argv[3]
split_rexe = sys.argv[4]
# split_rexe = 'chr[1-7]C'
gene_prex = sys.argv[5]
output_name = sys.argv[6]

fasta_dict = {}

output_bed = open(output_name + ".bed", "w")


with open(input_bed, "r") as bed_file:
    for line in bed_file:
        line = line.rstrip("\n")
        line_split = line.split("\t")
        chr_name = line_split[0].strip()
        gene_name = line_split[1]
        start = line_split[2]
        end = line_split[3]
        # 正则匹配
        # print(split_rexe)
        if re.search(split_rexe, chr_name) is None:
            continue
        else:
            # print(chr_name)
            tag = re.search(split_rexe, chr_name)[0]
            name = gene_name
            new_name = gene_prex + name
            fasta_dict[name] = new_name
            output_bed.write(chr_name + "\t" + new_name +
                             "\t" + start + "\t" + end + "\n")

output_bed.close()

cds_dict = {}
pep_dict = {}
with open(input_cds, "r") as fasta_file:
    for line in fasta_file:
        line = line.rstrip("\n")
        if line.startswith(">"):
            name = line[1:]
            cds_dict[name] = ""
        else:
            line = line + "\n"
            cds_dict[name] += line


with open(input_pep, "r") as fasta_file:
    for line in fasta_file:
        line = line.rstrip("\n")
        if line.startswith(">"):
            name = line[1:]
            pep_dict[name] = ""
        else:
            line = line + "\n"
            pep_dict[name] += line

output_cds = open(output_name + ".cds", "w")
output_pep = open(output_name + ".pep", "w")

for key in fasta_dict.keys():
    try:
        pep_dict[key]
        output_cds.write(">" + fasta_dict[key] + "\n")
        output_cds.write(cds_dict[key])
        output_pep.write(">" + fasta_dict[key] + "\n")
        output_pep.write(pep_dict[key])
    except KeyError:
        continue


output_cds.close()
output_pep.close()
