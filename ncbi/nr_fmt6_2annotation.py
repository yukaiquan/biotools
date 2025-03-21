#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import time
import subprocess
import re
import argparse


def argument_parser():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-i','--input', type=str, help='nr blast outfmt 6')
    parser.add_argument('-n','--nr_id', type=str, help='nr annotation file(plant_class_nrid.tsv)')
    parser.add_argument('-o','--output', type=str, help='output file')

    return parser.parse_args()


def main():
    args = argument_parser()
    input_file = args.input
    nr_id = args.nr_id
    output_file = args.output

    nr_mapping = {}
    with open(nr_id,'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#'):
                continue
            else:
                line_list = line.split()
                nr_id = line_list[0]
                # 第一列之后所有列
                nr_class = ' '.join(line_list[1:])
                if nr_id not in nr_mapping:
                    nr_mapping[nr_id] = nr_class
                else:
                    print('duplicate nr_id: ' + nr_id)


    flag_list = []
    out_list = []
    with open(input_file,'r') as f:
        for line in f:
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
                    if subject not in nr_mapping:
                        print('query not in nr_mapping: ' + subject)
                        class_name = "NA"
                    else:
                        class_name = nr_mapping[subject]
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
        sys.stderr.write("User interrupt me! ;-) See you!\n"
                           "%s" % time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
                           + "\n")
        sys.exit(0)