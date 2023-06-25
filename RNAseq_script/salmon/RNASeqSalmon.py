#!/usr/bin/env python3
import sys
import os
from subprocess import Popen, PIPE
# 多线程
from concurrent.futures import ProcessPoolExecutor
import shutil
import time


FASTP_OUT_DIR = '01_fastp'
SALMON_OUT_DIR = '02_salmon'
MATRIX_OUT_DIR = '03_matrix'
KALLISTO = '/public/home/acfaa2ssz7/soft/kallisto2matrix'


# 写出样本对应名称为matrix做准备
def writer_samples(samples: dict, genome_index: dict, tag: str) -> dict:
    sample_salmon_dict: dict = {}
    num: int = 0
    for key, value in samples.items():
        # 根据key.split('_')[1]写出不同的文件
        for key1, value1 in genome_index.items():
            if key.split('_')[-1] == key1:
                file_path = MATRIX_OUT_DIR + '/' + key1 + '_' + tag + '_samples.csv'
                # by01_1_sfs 改为 by01_1
                # 去掉最后的sfs
                key_list = key.split('_')
                key_list.pop()
                key = '_'.join(key_list)
                # 不存在则创建文件 存在则追加
                if not os.path.exists(file_path) or num == 0:
                    file_write = open(file_path, 'w')
                    file_write.write(value + ',' + key + '\n')
                    file_write.close()
                else:
                    file_write = open(file_path, 'a')
                    file_write.write(value + ',' + key + '\n')
                    file_write.close()
                sample_salmon_dict[file_path] = key1 + '_' + tag
                num += 1
    return sample_salmon_dict


# 判断文件夹和文件是否存在
def check_file(file_path: str) -> None:
    if not os.path.exists(file_path):
        print('Error: file or directory not exists ' + file_path)
        sys.exit(1)
    return None


# 判断文件夹是否存在，不存在则创建文件夹
def check_dir(dir_path: str) -> None:
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
    return None


# 并行数量计算
def cacu_threads(threads: int, intheads: int) -> int:
    if threads <= intheads:
        parallel = 1
        intheads = threads
    elif threads > intheads and threads % intheads == 0:
        parallel = threads // intheads
    else:
        parallel = threads // intheads + 1
    return parallel


def run_cmd(cmd: str) -> None:
    print(cmd)
    output = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    stdout, stderr = output.communicate()
    if output.returncode != 0:
        print(stderr)
        sys.exit(1)
    output.wait()
    return None


# fastp
def fastp(fastq_1: str, fastq_2: str, threads: int) -> list:
    cmd = 'fastp -i ' + fastq_1 + ' -I ' + fastq_2 + ' -o ' + FASTP_OUT_DIR + '/' + fastq_1.split('/')[-1].split('.')[0] + '_clean.fastq.gz -O ' + FASTP_OUT_DIR + '/' + fastq_2.split('/')[-1].split(
        '.')[0] + '_clean.fastq.gz -h ' + FASTP_OUT_DIR + '/' + fastq_1.split('/')[-1].split('.')[0] + '.html -j ' + FASTP_OUT_DIR + '/' + fastq_1.split('/')[-1].split('.')[0] + '.json' + ' -w ' + str(threads)
    print(cmd)
    output = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    stdout, stderr = output.communicate()
    if output.returncode != 0:
        print(stderr)
        sys.exit(1)
    output.wait()
    return [FASTP_OUT_DIR + '/' + fastq_1.split('/')[-1].split('.')[0] + '_clean.fastq.gz', FASTP_OUT_DIR + '/' + fastq_2.split('/')[-1].split('.')[0] + '_clean.fastq.gz']

# salmon


def salmon(genome_index: str, genome_gff: str, fastq_1: str, fastq_2: str, threads: int, tag: str) -> list:
    output_dir = SALMON_OUT_DIR + '/' + \
        fastq_1.split('/')[-1].split('.')[0].split('_')[0] + '_' + tag
    cmd = 'salmon quant -i ' + genome_index + ' -g ' + genome_gff + ' --gcBias -l A --numBootstraps 100 -1 ' + fastq_1 + ' -2 ' + fastq_2 + \
        ' -p ' + str(threads) + ' -o ' + output_dir
    return [cmd, output_dir + '/quant.sf', output_dir + '/quant.genes.sf']

# 合并salmon结果
# kallisto2matrix -i sfs_samples.csv -o sfs_tx -t salmon


def kallisto2matrix(samples: str, output_name: str) -> None:
    cmd = KALLISTO + ' -i ' + samples + ' -o ' + output_name + ' -t salmon'
    print(cmd)
    output = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    stdout, stderr = output.communicate()
    if output.returncode != 0:
        print(stderr)
        sys.exit(1)
    output.wait()
    return None


start_time = time.time()
# input salmon index
input_genome_index: str = sys.argv[1]
# 文件内容
# /public/home/acfaa2ssz7/genome/sfs/salmon/ACD_sativa_sfs_cds.fasta.index  /public/home/acfaa2ssz7/genome/sfs/salmon/ACD_sativa_sfs.gff    sfs
# input fastq file
input_fastq_file: str = sys.argv[2]
# input threads
threads: int = int(sys.argv[3])
# 文件内容
# SRR2201923_1.fastq.gz SRR2201923_2.fastq.gz   by01_1
# SRR2201924_1.fastq.gz SRR2201924_2.fastq.gz   by01_2
# SRR2201925_1.fastq.gz SRR2201925_2.fastq.gz   by02_1


genome_index: dict = {}

with open(input_genome_index, 'r') as f:
    for line in f:
        line = line.strip()
        if line.startswith('#'):
            continue
        else:
            # genome_index[line.split()[1]] = line.split()[0]
            if line.split()[1] in genome_index.keys():
                print('Error: duplicated key in genome index file')
                sys.exit(1)
            else:
                try:
                    check_file(line.split()[0])
                except:
                    print('Error: genome index file not exists')
                    sys.exit(1)
                genome_index[line.split()[2]] = [line.split()
                                                 [0], line.split()[1]]


check_dir(FASTP_OUT_DIR)
fastq_file_dict: dict = {}

with open(input_fastq_file, 'r') as f:
    for line in f:
        line = line.strip()
        if line.startswith('#'):
            continue
        else:
            if line.split()[2] in fastq_file_dict.keys():
                print('Error: duplicated key in fastq file')
                sys.exit(1)
            else:
                try:
                    check_file(line.split()[0])
                    check_file(line.split()[1])
                except:
                    print('Error: fastq file not exists')
                    sys.exit(1)
                fastq_file_dict[line.split()[2]] = line.split()[
                    0] + ' ' + line.split()[1]


# 运行fastp
fastp_threads: int = 16
fastp_parallel: int = cacu_threads(threads, fastp_threads)
fastp_out_dict: dict = {}
with ProcessPoolExecutor(max_workers=fastp_parallel) as executor:
    for key, value in fastq_file_dict.items():
        fastq_1 = value.split()[0]
        fastq_2 = value.split()[1]
        fastp_out_dict[key] = executor.submit(
            fastp, fastq_1, fastq_2, fastp_threads)

# 运行salmon
salmon_threads: int = 16
check_dir(SALMON_OUT_DIR)
# 先生成salmon命令
salmon_cmd_list: list = []
salmon_gene_out_dict: dict = {}
salmon_transcript_out_dict: dict = {}
for key, value in genome_index.items():
    for key2, value2 in fastp_out_dict.items():
        fastq_1 = value2.result()[0]
        fastq_2 = value2.result()[1]
        # salmon_cmd_list.append(
        #     salmon(value[0], value[1], fastq_1, fastq_2, salmon_threads, key)[0])
        # salmon_gene_out_dict[key2 + '_' + key] = salmon(
        #     value[0], value[1], fastq_1, fastq_2, salmon_threads, key)[1]
        salmon_cmd, salmon_gene_out, salmon_transcript_out = salmon(
            value[0], value[1], fastq_1, fastq_2, salmon_threads, key)
        salmon_cmd_list.append(salmon_cmd)
        salmon_gene_out_dict[key2 + '_' + key] = salmon_gene_out
        salmon_transcript_out_dict[key2 + '_' + key] = salmon_transcript_out

# 再运行salmon
salmon_parallel: int = cacu_threads(threads, salmon_threads)
with ProcessPoolExecutor(max_workers=salmon_parallel) as executor:
    for cmd in salmon_cmd_list:
        executor.submit(run_cmd, cmd)

# 检查文件夹是否存在，不存在则创建
check_dir(MATRIX_OUT_DIR)
# sfs_tx_samples.csv
# 02_salmon_out/SRR22937062_sfs/quant.sf,by01_1
# 02_salmon_out/SRR22937061_sfs/quant.sf,by01_2
# 02_salmon_out/SRR22937058_sfs/quant.sf,by01_3

# print(salmon_gene_out_dict)
# print(salmon_transcript_out_dict)

sample_salmon_gene_dict: dict = writer_samples(
    salmon_gene_out_dict, genome_index, 'gene')
sample_salmon_transcript_dict: dict = writer_samples(
    salmon_transcript_out_dict, genome_index, 'transcript')

# print(sample_salmon_gene_dict)
# print(sample_salmon_transcript_dict)

# 合并gene和transcript
salmon_samples_list: list = []
for key, value in sample_salmon_gene_dict.items():
    salmon_samples_list.append([key, value])
for key, value in sample_salmon_transcript_dict.items():
    salmon_samples_list.append([key, value])

# 运行kallisto2matrix
kallisto2matrix_parallel: int = cacu_threads(threads, 1)
with ProcessPoolExecutor(max_workers=kallisto2matrix_parallel) as executor:
    for sample in salmon_samples_list:
        executor.submit(kallisto2matrix, sample[0], sample[1])

# 将结果文件移动到指定文件夹
# 提取salmon_samples_list第二列去重
salmon_samples_list_2: list = []
for sample in salmon_samples_list:
    salmon_samples_list_2.append(sample[1])
salmon_samples_list_2 = list(set(salmon_samples_list_2))
# 将列表中以此为文件前缀的文件移动到MATRIX_OUT_DIR
for sample in salmon_samples_list_2:
    if os.path.exists(sample + '_count_matrix.txt'):
        shutil.move(sample + '_count_matrix.txt', MATRIX_OUT_DIR)
    else:
        print('Error: count matrix file not exists')
    if os.path.exists(sample + '_tpm_matrix.txt'):
        shutil.move(sample + '_tpm_matrix.txt', MATRIX_OUT_DIR)
    else:
        print('Error: tpm matrix file not exists')
    if os.path.exists(sample + '_fpkm_matrix.txt'):
        shutil.move(sample + '_fpkm_matrix.txt', MATRIX_OUT_DIR)
    else:
        print('Error: fpkm matrix file not exists')


print('time used: ' + str(time.time() - start_time))
print('##############################All done!########################################')
