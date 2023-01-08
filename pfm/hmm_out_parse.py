#!/use/bin/env pyhton3

import numpy as np
import gzip
import time
from tqdm import tqdm


input_hmm_out_list = ['sfs.output.gz', 'SAU_ins.fasta.output.gz',
                      'SAU_lon.fasta.output.gz', 'Sang_lon.fasta.output.gz',
                      'strigosa.fasta.output.gz', 'atlantica.fasta.output.gz',
                      'OTv2.fasta.output.gz', 'Sang.fasta.output.gz',
                      'eriantha.fasta.output.gz', 'Sang_ins.fasta.output.gz']


def hmm_out_parse(input_hmm_out):
    """
    :param input_hmm_out: hmmsearch的结果文件
    :return: hmmsearch的结果文件
    """
    hmmout_list = []
    # 读取压缩文件
    with gzip.open(input_hmm_out, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            else:
                # print(line.strip().split())
                line = line.strip().split()
                gene_name = line[0]
                query_name = line[3]
                accession = line[4]
                e_value = line[6]
                score = line[7]
                alig_start = line[17]
                alig_end = line[18]
                hmmout_list.append(
                    [gene_name, query_name, accession, e_value, score, alig_start, alig_end])

    # 将列表转换为数组
    hmmout_array = np.array(hmmout_list)
    # 按照第一列进行排序
    hmmout_array = hmmout_array[hmmout_array[:, 0].argsort()]
    # 根据e_value<=1e-5进行筛选
    hmmout_array = hmmout_array[hmmout_array[:, 3].astype(float) <= 1e-5]

    return hmmout_array


# 合并所有的hmmsearch结果
hmmout_array = np.array([])
for input_hmm_out in tqdm(input_hmm_out_list, desc='hmm_out_parse'):
    hmmout_array = np.vstack((hmmout_array, hmm_out_parse(input_hmm_out))) if hmmout_array.size else hmm_out_parse(
        input_hmm_out)

# 添加目前时间
now_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
# 最后添加两列时间
hmmout_array = np.column_stack(
    (hmmout_array, np.full((hmmout_array.shape[0], 2), now_time)))
# 在第一列添加一列序号
hmmout_array = np.column_stack(
    (np.arange(1, hmmout_array.shape[0] + 1), hmmout_array))
# 写出文件
np.savetxt('hmm_out_parse.txt', hmmout_array, fmt='%s', delimiter='\t')
