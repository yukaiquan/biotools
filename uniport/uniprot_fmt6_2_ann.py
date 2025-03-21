#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import os
import time
import subprocess
import re
import argparse
import tqdm


def argument_parser():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-i','--input', type=str, help='uniprot blast outfmt 6')
    parser.add_argument('-u','--uniprot', type=str, help='uniprot file')
    parser.add_argument('-o','--output', type=str, help='output file')

    return parser.parse_args()

def main():
    args = argument_parser()
    input_file = args.input
    uniprot = args.uniprot
    output_file = args.output

    uniprot_dict = {}
    with open(uniprot, 'r') as f:
        for line in tqdm.tqdm(f, desc='reading uniprot'):
            if line.startswith('#'):
                continue
            lines = line.strip().split('\t')
            id = lines[0]
            if id not in uniprot_dict:
                # print(lines)
                if len(lines) == 2:
                    uniprot_dict[id] = lines[1]
                else:
                    uniprot_dict[id] = lines[2]
            else:
                print('duplicate id:', id)
        flag_list = []
    out_list = []
    with open(input_file,'r') as f:
        for line in tqdm.tqdm(f, desc='reading blast'):
            line = line.strip()
            if line.startswith('#'):
                continue
            else:
                line_list = line.split('\t')
                query = line_list[0]
                if query not in flag_list:
                    subject = line_list[1]
                    identical = line_list[2]
                    evalue = line_list[10]
                    score = line_list[11]
                    class_name = ""
                    if subject not in uniprot_dict:
                        print('query not in uniprot: ' + subject)
                        class_name = "NA"
                    else:
                        class_name = uniprot_dict[subject]
                    out_list.append([query,subject,identical,evalue,score,class_name])
                else:
                    continue
                flag_list.append(query)
    with open(output_file,'w') as f:
        f.write('query\tsubject\tidentical\tevalue\tscore\tannotation\n')
        for line in out_list:
            f.write('\t'.join(line) + '\n')

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me! ;-) See you!\n")
        sys.exit(0)