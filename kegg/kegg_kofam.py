#! /usr/bin/env python3
import gzip
# import pandas as pd


blastp_file_list: list = []
score: float = 100
evalue: float = 1e-5
with gzip.open('ACD.SFS.pep.fasta.querry2KO.gz', 'rb') as f:
    for line in f:
        line = line.decode('utf-8')
        if line.startswith('#'):
            continue
        else:
            line = line.strip().replace('*', '').split('\t')
            print(line)
            if float(line[4]) >= score and float(line[5]) <= evalue:
                blastp_file_list.append(
                    [line[1], line[2], line[3], line[4], line[5], line[6]])
    # convert list to ndarray

input_kegg[0].str.split(' ', expand=True)
