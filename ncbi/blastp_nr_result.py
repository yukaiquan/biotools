'''
yukaiquan 1962568272@qq.com
2022-09-14
A.satnudSFSUN01G001420	XP_047045654.1	68.9	315	55	7	54	334	5	310	1.8e-104	385.6
identity 80 align_len 100 evalue 1e-5
'''
from tqdm import tqdm  # progress bar
import numpy as np
import pandas as pd
import os
import sys

def main(blastp_out:str, nr_id_function_filename:str,identity:float,align_len:float,evalue:float,out_put:str):
    blastp_file_list: list = []

    with open(blastp_out, 'r') as f:
        for line in tqdm(f, desc='Read blastp result file and filter'):
            if line.startswith('#'):
                continue
            else:
                line = line.strip().split('\t')
                if float(line[2]) >= identity and float(line[3]) >= align_len and float(line[10]) <= evalue:
                    blastp_file_list.append(
                        [line[0], line[1], line[2], line[3], line[10], line[11]])
    # convert list to ndarray
    blastp_file_ndarray = np.array(blastp_file_list)
    # sort by identity and name
    blastp_file_ndarray = blastp_file_ndarray[np.lexsort(
        (blastp_file_ndarray[:, 0], blastp_file_ndarray[:, 5][::-1]))]
    # ES6 array to list by name
    vals, idx_start, count = np.unique(
        blastp_file_ndarray[:, 0], return_counts=True, return_index=True)
    # extract the first genename
    blastp_file_ndarray = blastp_file_ndarray[idx_start]

    nr_function = {}
    with open(nr_id_function_filename, "rb") as f:
        for fLine in tqdm(f, desc='Read nr id function file'):
            fLine = fLine.decode().strip().split()
            nr_function[fLine[0]] = ' '.join(fLine[1:])

    function_ndarry = np.array(list(map(lambda x: nr_function.get(x,None), list(blastp_file_ndarray[:, 1]))))
    blastp_file_ndarray = pd.DataFrame(blastp_file_ndarray)
    nr_function = pd.DataFrame(function_ndarry)
    blastp_file_ndarray = pd.concat([blastp_file_ndarray, nr_function], axis=1)
    blastp_file_ndarray.to_csv(out_put, sep='\t', index=False, header=False)

if __name__ == '__main__':
    try:
        blastp_out: str = sys.argv[1]
        nr_id_function_filename: str = sys.argv[2]
        identity: float = float(sys.argv[3]) # 80
        align_len: float = float(sys.argv[4]) # 100
        evalue: float = float(sys.argv[5])  # 1*10^-5
        out_put: str = sys.argv[6]
        print('blastp_2_ncbi:'+ blastp_out + 'nr_id_function_filename:xp1928.1 function:'+ nr_id_function_filename+ 'identity:'  + str(identity) +  'align_len:' + str(align_len) + 'evalue:' + str(evalue) + 'out_put:' + out_put)
        main(blastp_out, nr_id_function_filename,identity,align_len,evalue,out_put)
    except KeyboardInterrupt:
        print('Interrupted')
        os._exit(0)
 