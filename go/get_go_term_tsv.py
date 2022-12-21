#! /use/bin/env python3
'''
python get_go_term_sql.py input_file output_file
go-basic.obo turn to go_term.tsv
mysql namespace:
    biological_process 1
    molecular_function 2
    cellular_component 3
go_term.tsv format:
go_id	name	namespace	created_time	updated_time
GO:0000001	mitochondrion inheritance	1	2022-09-30 09:53:15	2022-09-30 09:53:15
GO:0000002	mitochondrial genome maintenance	1	2022-09-30 09:53:15	2022-09-30 09:53:15
GO:0000003	reproduction	1	2022-09-30 09:53:15	2022-09-30 09:53:15
GO:0000005	obsolete ribosomal chaperone activity	2	2022-09-30 09:53:15	2022-09-30 09:53:15
GO:0000006	high-affinity zinc transmembrane transporter activity	2	2022-09-30 09:53:15	2022-09-30 09:53:15
GO:0000007	low-affinity zinc ion transmembrane transporter activity	2	2022-09-30 09:53:15	2022-09-30 09:53:15
'''

import sys
import pandas as pd


biological_process: str = '1'
molecular_function: str = '2'
cellular_component: str = '3'

term_list: list = []
input_file = sys.argv[1]
output_file = sys.argv[2]
# input_file = "go-basic.obo"
# output_file = "go_term.tsv"
raw_file = open(input_file).read()

go_alt = {}
for go_term in raw_file.split("[Term]"):
    go_id = ''
    name = ''
    namespace = ''
    for line in go_term.split("\n"):
        if str(line).startswith('[Typedef]'):
            break
        if str(line).startswith("id") and "GO:" in line:
            go_id = line.rstrip().split(" ")[1]
        if str(line).startswith("name:"):
            name = line.rstrip().split(": ")[1]
        if str(line).startswith("namespace"):
            namespace = line.rstrip().split(" ")[1].strip()
        if str(line).startswith("alt_id"):
            # 添加alt_id 它的值和go_id一样
            alt_id = line.rstrip().split(" ")[1]
            go_alt[alt_id] = go_id
    term = go_id + '\t' + name + '\t' + namespace
    if '' != go_id:
        term_list.append(term)
        for key in go_alt.keys():
            if go_id == go_alt[key]:
                term = key + '\t' + name + '\t' + namespace
                term_list.append(term)

term_df = pd.DataFrame(term_list)
term_df = term_df[0].str.split('\t', expand=True)
term_df.columns = ['go_id', 'name', 'namespace']

term_df['name'] = term_df['name'].str.replace('\"', '')
term_df['namespace'] = term_df['namespace'].str.replace(
    'biological_process', biological_process)
term_df['namespace'] = term_df['namespace'].str.replace(
    'molecular_function', molecular_function)
term_df['namespace'] = term_df['namespace'].str.replace(
    'cellular_component', cellular_component)

created_time = pd.to_datetime('today').strftime('%Y-%m-%d %H:%M:%S')
term_df['created_time'] = created_time
term_df['updated_time'] = created_time
# 根据go_id去排序
term_df = term_df.sort_values(by='go_id')

term_df.to_csv(output_file, index=False, header=True, sep='\t')
