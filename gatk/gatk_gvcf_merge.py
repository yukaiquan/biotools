import os
import sys
from subprocess import Popen, PIPE
import argparse
import time
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor

BCFTOOLS_PATH = '/usr/local/bin/bcftools'


def bcftools_index(vcf_file: str) -> bool:
    """
    bcftools must be installed 1.3.1
    :param vcf_file:
    :return:
    """
    if not file_exists(vcf_file):
        print('vcf file not exists')
        return False
    cmd = '{BCFTOOLS_PATH} index -f -t {vcf_file}'.format(
        BCFTOOLS_PATH=BCFTOOLS_PATH, vcf_file=vcf_file)
    pop = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    pop.wait()
    if pop.returncode != 0:
        print(pop.stderr.read().decode())
        return False
    return True


def file_exists(file_name: str) -> bool:
    """
    判断文件是否存在当前目录
    :param file_path:
    :return:
    """
    if os.path.exists(file_name):
        return True
    return


def gatk_CombineGVCFs(ref_file: str, gvcf_list: list, out_file: str) -> bool:
    '''
    :param gvcf_list: gvcf文件列表
    :param out_file: 输出文件
    :return:
    '''
    if not isinstance(gvcf_list, list):
        print('gvcf_list must be a list')
        return False
    if file_exists(out_file):
        print('out file exists')
        return False
    cmd = 'gatk CombineGVCFs -OVI False -R {ref_file} -O {out_file} -V {gvcf_list}'
    cmd = cmd.format(ref_file=ref_file, out_file=out_file,
                     gvcf_list=' -V '.join(gvcf_list))
    pop = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    pop.wait()
    return True


def gatk_GenotypeGVCFs(ref_file: str, gvcf_file: str, out_file: str) -> bool:
    '''
    :param gvcf_file: gvcf文件
    :param out_file: 输出文件
    :return:
    '''
    if not file_exists(gvcf_file):
        print('gvcf file not exists')
        return False
    if file_exists(ref_file):
        print('referrence file exists')
        return False
    cmd = 'gatk GenotypeGVCFs -OVI False -R {ref_file} -O {out_file} -V {gvcf_file}'
    cmd = cmd.format(ref_file=ref_file, out_file=out_file,
                     gvcf_file=gvcf_file)
    pop = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    pop.wait()
    return True


def main():
    # 进度条
    tqdm.monitor_interval = 0
    start_time = time.time()
    # 获取参数
    parser = argparse.ArgumentParser(description='gatk gvcf merge to vcf')
    parser.add_argument('-r', '--ref', help='referrence file', required=True)
    parser.add_argument('-g', '--gvcf', help='gvcf file list', required=True)
    parser.add_argument(
        '-o', '--out', help='out file name prefix', required=True)
    parser.add_argument(
        '-t', '--thread', help='thread number', default=1, type=int)
    args = parser.parse_args()
    ref_file = args.ref
    gvcf_list_file = args.gvcf
    out_file = args.out
    thread_num = args.thread
    # 读取gvcf文件列表
    gvcf_list = []
    if not file_exists(gvcf_list_file):
        print('gvcf list file not exists')
        sys.exit(1)
    with open(gvcf_list_file, 'r') as f:
        for line in f:
            gvcf_list.append(line.strip())
    # 判断索引文件是否存在
    with ProcessPoolExecutor(max_workers=thread_num) as executor:
        for gvcf_file in gvcf_list:
            executor.submit(bcftools_index, gvcf_file)
    tqdm.write('index file check done')
    # 合并gvcf文件
    if not gatk_CombineGVCFs(ref_file, gvcf_list, out_file + '.g.vcf.gz'):
        print('gvcf merge failed')
        sys.exit(1)
    tqdm.write('gvcf merge done')
    if not bcftools_index(out_file + '.g.vcf.gz'):
        print('index file create failed')
        sys.exit(1)
    tqdm.write('index file create done')
    # 生成vcf文件
    if not gatk_GenotypeGVCFs(ref_file, out_file + '.g.vcf.gz', out_file + '.vcf.gz'):
        print('gvcf merge failed')
        sys.exit(1)
    tqdm.write('vcf file create done')
    if not bcftools_index(out_file + '.vcf.gz'):
        print('index file create failed')
        sys.exit(1)
    tqdm.write('index file create done')
    end_time = time.time()
    print('time used: {}s'.format(end_time - start_time))


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me! ;-) Bye!")
