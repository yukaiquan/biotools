#!/use/bin/env python3

import sys

input_list: str = sys.argv[1]
output_file: str = sys.argv[2]

list_file: list = []

with open(input_list, 'r') as f:
    line = f.readlines()
    for l in line:
        if l not in list_file:
            list_file.append(l)
        else:
            print(l)

with open(output_file, 'w') as f:
    for l in list_file:
        f.write(l)
