'''
yukaiquan 1962568272@qq.com
2022-09-14
A.satnudSFSUN01G001420	XP_047045654.1	68.9	315	55	7	54	334	5	310	1.8e-104	385.6
identity 60 align_len 100 evalue 1e-5
'''
from tqdm import tqdm  # progress bar
import numpy as np
import pandas as pd
import os
import sys

def nr_id_function(nr_id_function_filename:str)->dict:
    nr_function = {}
    with open(nr_id_function_filename, "rb") as f:
        for fLine in tqdm(f, desc='Read nr id function file'):
            fLine = fLine.decode().strip().split(',')
            nr_function[fLine[0]] = ' '.join(fLine[1:-2])
    return nr_function

def writer_nrout(nr_function:dict, blastp_out_file:str, nr_output:str,identity:float, align_len:float, evalue:float)->None:
    blastp_file_list: list = []
    with open(blastp_out_file, 'r') as f:
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



    # memory is to large
    function_ndarry = np.array(list(map(lambda x: nr_function.get(x,None), list(blastp_file_ndarray[:, 1]))))
    
    
    blastp_file_ndarray = pd.DataFrame(blastp_file_ndarray)
    nr_functions = pd.DataFrame(function_ndarry)
    print("start merge")
    # blastp_file_nd = np.concatenate((blastp_file_ndarray, function_ndarry), axis=1)
    
    blastp_file_ndarray = pd.concat([blastp_file_ndarray, nr_functions], axis=1)
    # np.savetxt(nr_output, blastp_file_nd, fmt='%s', delimiter='\t')
    blastp_file_ndarray.to_csv(nr_output, sep=',', index=False, header=False)
    del blastp_file_ndarray
    del nr_functions
    del function_ndarry
    del blastp_file_list

def main( nr_id_function_filename:str,identity:float,align_len:float,evalue:float,blastp_outs:str)->None:
    nr_function = nr_id_function(nr_id_function_filename)
    for blastp_out in blastp_outs.split(','):
        nr_output = blastp_out + '.nrout'
        writer_nrout(nr_function, blastp_out, nr_output,identity,align_len,evalue)
    
 




if __name__ == '__main__':
    try:
        nr_id_function_filename: str = sys.argv[1]
        identity: float = float(sys.argv[2]) # 80
        align_len: float = float(sys.argv[3]) # 100
        evalue: float = float(sys.argv[4])  # 1*10^-5
        blastp_outs: str = sys.argv[5]
        print('blastp_2_ncbi:'+ blastp_outs + 'function:'+ nr_id_function_filename+ 'identity:'  + str(identity) +  'align_len:' + str(align_len) + 'evalue:' + str(evalue) )
        main(nr_id_function_filename,identity,align_len,evalue,blastp_outs)
    except KeyboardInterrupt:
        print('Interrupted')
        os._exit(0)
 