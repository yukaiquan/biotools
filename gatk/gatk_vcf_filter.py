#!/usr/bin/env python

import os
import sys
from subprocess import Popen, PIPE
import argparse
import time
# 多线程
from concurrent.futures import ProcessPoolExecutor

BCFTOOLS_PATH = '/public/home/pengyuanying/.conda/envs/gatk42/bin/bcftools'


def file_exists(file_name:str) -> bool:
    """
    判断文件是否存在当前目录
    :param file_path:
    :return:
    """
    if os.path.exists(file_name):
        return True
    return False


def bcftools_index(vcf_file) -> bool:
    """
    bcftools must be installed 1.3.1
    :param vcf_file:
    :return:
    """
    cmd = '{BCFTOOLS_PATH} index -f -t {vcf_file}'.format(
        BCFTOOLS_PATH=BCFTOOLS_PATH, vcf_file=vcf_file)
    pop = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    pop.wait()
    if pop.returncode != 0:
        print(pop.stderr.read().decode())
        return False
    return True

def bcftools_index_csi(vcf_file) -> bool:
    """
    bcftools must be installed 1.3.1
    :param vcf_file:
    :return:
    """
    cmd = '{BCFTOOLS_PATH} index -c {vcf_file}'.format(
        BCFTOOLS_PATH=BCFTOOLS_PATH, vcf_file=vcf_file)
    pop = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    pop.wait()
    if pop.returncode != 0:
        print(pop.stderr.read().decode())
        return False
    return True

def gatk_SelectVariants(vcf_file, selet_type, out_file) -> bool:
    '''
    :param vcf_file: vcf文件
    :param selet_type: 选择类型
    :param out_file: 输出文件
    :return:
    '''
    if not file_exists(vcf_file):
        print('vcf file not exists')
        return False
    if file_exists(out_file):
        print('out file exists')
        return False
    cmd = 'gatk SelectVariants -OVI False -V {vcf_file} -select-type {selet_type} -O {out_file}'
    cmd = cmd.format(vcf_file=vcf_file,
                     selet_type=selet_type, out_file=out_file)
    pop = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    pop.wait()
    if pop.returncode != 0:
        print(pop.stderr.read().decode())
        return False
    return True


def gatk_VariantFiltration(vcf_file, out_file, snp=True) -> bool:
    '''
    :param vcf_file: vcf文件
    :param out_file: 输出文件
    :param snp: 是否是snp 默认是
    :return:
    '''
    if not file_exists(vcf_file):
        print('vcf file not exists')
        return False
    if file_exists(out_file):
        print('out file exists')
        return False
    if snp:
        cmd = 'gatk VariantFiltration -OVI False -V {vcf_file}  -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "SOR > 3.0" --filter-name "SOR3" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O {out_file}'
    else:
        cmd = 'gatk VariantFiltration -OVI False -V {vcf_file} -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "FS > 200.0" --filter-name "FS200" -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" -O {out_file}'
    cmd = cmd.format(vcf_file=vcf_file, out_file=out_file)
    # 忽略警告
    pop = Popen(cmd, shell=True)
    pop.wait()
    # print(pop.stderr.read().decode())
    return True


def gatk_VariantsToTable(vcf_file, out_file, snp=True) -> bool:
    '''
    :param vcf_file: vcf文件
    :param out_file: 输出文件
    :param snp: 是否是snp 默认是
    :return:
    '''
    if not file_exists(vcf_file):
        print('vcf file not exists')
        return False
    if file_exists(out_file):
        print('out file exists')
        return False
    if snp:
        cmd = 'gatk VariantsToTable -V {vcf_file} -F CHROM -F POS -F REF -F ALT -F FILTER -F DP -F QUAL -F QD -F SOR -F FS -F MQ -O {out_file}'
    else:
        cmd = 'gatk VariantsToTable -V {vcf_file} -F CHROM -F POS -F REF -F ALT -F FILTER -F DP -F QUAL -F QD -F SOR -F FS -F MQ -O {out_file}'
    cmd = cmd.format(vcf_file=vcf_file, out_file=out_file)
    pop = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    pop.wait()
    if pop.returncode != 0:
        print(pop.stderr.read().decode())
        return False
    return True


def main():
    start_time = time.time()
    parser = argparse.ArgumentParser(description='gatk vcf filter')
    parser.add_argument('-i', '--input', help='input vcf file', required=True)
    parser.add_argument(
        '-o', '--output', help='output file name prefix', required=True)
    args = parser.parse_args()
    vcf_file = args.input
    out_file = args.output
    # 判断vcf文件和索引文件是否存在
    print('check vcf file and index file')
    if not file_exists(vcf_file):
        print('vcf file not exists')
        sys.exit(1)
    if not file_exists(vcf_file + '.tbi'):
        if not bcftools_index(vcf_file):
            print('bcftools index error')
            sys.exit(1)
    # 判断输出文件是否存在
    # 多线程同时运行
    # 创建两个线程
    file_list = [['SNP', out_file + '.snp.vcf.gz'],
                 ['INDEL', out_file + '.indel.vcf.gz']]
    filter_list = [[out_file + '.snp.vcf.gz', out_file + '.snp.filter.vcf.gz', True],
                   [out_file + '.indel.vcf.gz', out_file + '.indel.filter.vcf.gz', False]]
    table_list = [[out_file + '.snp.filter.vcf.gz', out_file + '.snp.filter.txt', True],
                  [out_file + '.indel.filter.vcf.gz', out_file + '.indel.filter.txt', False]]
    print('start gatk SelectVariants')
    with ProcessPoolExecutor(max_workers=2) as executor:
        for snp, file in file_list:
            if not file_exists(file):
                executor.submit(gatk_SelectVariants, vcf_file, snp, file)
    # 创建两个索引
    print('start bcftools index')
    if not file_exists(out_file + '.snp.vcf.gz.tbi') or not file_exists(out_file + '.snp.vcf.gz.csi'):
        with ProcessPoolExecutor(max_workers=2) as executor:
            for snp, file in file_list:
                if not file_exists(file + '.tbi') or not file_exists(file + '.csi'):
                    executor.submit(bcftools_index, file)
    # 过滤snp
    print('start gatk VariantFiltration')
    with ProcessPoolExecutor(max_workers=2) as executor:
        for infile, outfile, snp in filter_list:
            if not file_exists(outfile):
                executor.submit(gatk_VariantFiltration, infile, outfile, snp)
    # 创建两个索引
    print('start bcftools index')
    with ProcessPoolExecutor(max_workers=2) as executor:
        for infile, outfile, snp in filter_list:
            if not file_exists(outfile + '.tbi') or not file_exists(outfile + '.csi'):
                executor.submit(bcftools_index, outfile)
    # 导出表格
    print('start gatk VariantsToTable')
    with ProcessPoolExecutor(max_workers=2) as executor:
        for infile, outfile, snp in table_list:
            if not file_exists(outfile):
                executor.submit(gatk_VariantsToTable, infile, outfile, snp)
    end_time = time.time()
    print('time used: {}s'.format(end_time - start_time))


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print('interrupted by user, exiting ...')
        sys.exit(0)
