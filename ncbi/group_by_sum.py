#! /usr/bin/env python3

import pandas as pd
import sys

input_file = sys.argv[1]
input_file = "Avena_oat_all_ann_detail.tsv"
group_by = int(sys.argv[2])
sum_by = int(sys.argv[3])
output_file = sys.argv[4]

data_frame = pd.read_csv(input_file, sep="\t", header=None)
data_frame = data_frame.groupby([1]).agg({0: "sum"})
data_frame = data_frame.sort_values(by=[0], ascending=False).reset_index()
data_frame[3] = data_frame[1]
del data_frame[1]
data_frame.to_csv(output_file, sep="\t", header=None, index=None)
