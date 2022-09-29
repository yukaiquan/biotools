#!/use/bin/env python3

import pandas as pd
import sys
'''
python ko_list_2_sql.py ko_list.txt ko_list.csv
'''

# Read in the KO list
ko_list_input = sys.argv[1]
output_file = sys.argv[2]
# ko_list_input = r'D:\BDdatabase\KEGG\database\ko_list'
# output_file = r'D:\BDdatabase\KEGG\database\ko_list.csv'
ko_list = pd.read_csv(ko_list_input, sep='\t', header=0)

# add time column
ko_list['created_time'] = '2022-09-29 00:00:00'
ko_list['updated_time'] = '2022-09-29 00:00:00'

# write to csv
ko_list.to_csv(output_file, index=False, header=True)
