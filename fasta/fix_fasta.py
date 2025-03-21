#!/usr/bin/env python

import sys
import argparse

input_fasta = sys.argv[1]
output_fasta = sys.argv[2]

res_dict = {}

name = ""
with open(input_fasta, "r") as f:
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            name = line.split()[0]
            if name in res_dict:
                print(f"Duplicate header: {line}")
                sys.exit(1)
            res_dict[name] = ""
        else:
            line = line.replace(" ", "").replace("-", "").replace(".", "").replace("*", "").replace("_", "")
            if name in res_dict:
                res_dict[name] += line
            else:
                print(f"Header not found: {line}")
                sys.exit(1)

out_open = open(output_fasta, "w")

for name in res_dict.keys():
    out_open.write(f"{name}\n")
    # out_open.write(f"{res_dict[name]}\n")
    # 60个字符一行
    out_open.write("\n".join([res_dict[name][i:i+60] for i in range(0, len(res_dict[name]), 60)]))
    out_open.write("\n")

out_open.close()

        
print("Fasta file fixed successfully")