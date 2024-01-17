#!/usr/bin/env python3
import sys
import os
import time
import glob
from subprocess import Popen, PIPE
from concurrent.futures import ProcessPoolExecutor
import argparse

# python centromere_finding.py -g ~/genome/OT3098v2/ACD_sativa_OT3098v2_genome.fasta -i CRR515329_f1.fq.gz -I CRR515329_r2.fq.gz -c CRR515328_f1.fq.gz -C CRR515329_r2.fq.gz -t 16 -f otv2


FASTP = "01_fastp"
BAM = "02_bam"
MACS = "03_macs2"
CHIPLOT = "/public/home/acfaa2ssz7/soft/chip_plot"


def check_file(file_path: str) -> bool:
    '''
    检查文件是否存在
    '''
    return os.path.exists(file_path)


def mk_dir(dir_path: str) -> None:
    '''
    创建文件夹
    '''
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)


def cacu_threads(threads: int, intheads: int) -> int:
    '''
    计算并行数
    threads: 总线程数
    intheads: 每个并行任务的线程数
    '''
    if threads <= intheads:
        parallel = 1
        intheads = threads
    elif threads > intheads and threads % intheads == 0:
        parallel = threads // intheads
    else:
        parallel = threads // intheads + 1
    return parallel


def run_fastp(fastq1, fastq2, outname1, outname2, name, threads):
    # 测试cmd 睡眠10s
    # cmd = 'echo "hello world"\nsleep 2'
    # fastp最大线程数为16
    # if outname1 != outname2:
    #     print('Error: outname1 != outname2')
    #     sys.exit(1)
    cmd = f'fastp -i {fastq1} -I {fastq2} -o {FASTP}/{outname1} -O {FASTP}/{outname2} -h {FASTP}/{name}.html -j {FASTP}/{name}.json -w {threads}'
    p = Popen(cmd, shell=True)
    return_code = p.wait()  # 等待子进程结束，并返回状态码；
    return return_code


def run_bwa(ref: str, fastq1: str, fastq2: str, flag: str, threads: int) -> int:
    sample_title = '\'@RG\\tID:' + flag + '\\tSM:' + \
        flag + '\\tLB:' + flag + '\\tPL:Illumina\''
    cmd = f'bwa mem {ref} {FASTP + "/" +fastq1} {FASTP + "/" + fastq2} -t {threads} -R {sample_title} | samtools sort -@ {threads} -o {BAM +"/" +flag +".bam"}  && samtools view -@ {threads} -h -F 20 -o {BAM + "/" + flag + "_F20.bam"} {BAM +"/" +flag +".bam"}'
    p = Popen(cmd, shell=True)
    return_code = p.wait()
    print(cmd)
    if return_code != 0:
        print(f'Error: {cmd}')
        sys.exit(1)
    return return_code


def run_rmdup(bam: str, outbam: str) -> int:
    cmd = f'gatk --java-options -Xmx64G MarkDuplicates -I {BAM + "/" + bam} -O {BAM + "/" + outbam} -M {BAM + "/" + outbam + ".metrics.txt"} --ASSUME_SORTED true'
    p = Popen(cmd, shell=True)
    return_code = p.wait()
    print(cmd)
    if return_code != 0:
        print(f'Error: {cmd}')
        sys.exit(1)
    return return_code


def run_macs2(ck_bam: str, treat_bam: str, flag: str, genome_size: int) -> int:
    cmd = f'macs2 callpeak --nomodel --extsize 200 -f BAMPE --keep-dup all -c {BAM+"/" + ck_bam} -t {BAM+"/" + treat_bam} -n {flag} -g {genome_size} --outdir {MACS} --bdg -q 0.05 && macs2 bdgcmp -t {MACS + "/" + flag + "_treat_pileup.bdg" } -c {MACS + "/" + flag + "_control_lambda.bdg" } -o {MACS + "/" + flag + "_chip.bdg" } -m FE && sort -k1,1 -k2,2n {MACS + "/" + flag + "_chip.bdg" } | grep -i -v "chrun" > {MACS + "/" + flag + "_chip_sorted.bdg" }'
    p = Popen(cmd, shell=True)
    return_code = p.wait()
    print(cmd)
    if return_code != 0:
        print(f'Error: {cmd}')
        sys.exit(1)
    return return_code


def run_bedgraph(input: str, genomefai: str, output: str) -> int:
    cmd = f'bedGraphToBigWig {input} {genomefai} {output}'
    p = Popen(cmd, shell=True)
    return_code = p.wait()
    print(cmd)
    if return_code != 0:
        print(f'Error: {cmd}')
        sys.exit(1)
    return return_code


def chip_plot(flag: str, genomefai: str) -> int:
    # chip_plot -b sfs_chip_sorted.bdg -f ~/genome/sfs/bwaindex/ACD.sativa.fasta.fai -o sfs_chip_sorted.bin -s 100000 -w 1000000
    cmd = f'{CHIPLOT} -b {MACS + "/" + flag + "_chip_sorted.bdg"} -f {genomefai} -o {MACS + "/" + flag + "_chip_sorted.bin"} -s 100000 -w 1000000'
    p = Popen(cmd, shell=True)
    return_code = p.wait()
    print(cmd)
    if return_code != 0:
        print(f'Error: {cmd}')
        sys.exit(1)
    return return_code


def get_genome_size(genome: str) -> int:
    faidx_file = f'{genome}.fai'
    result = 0
    if not os.path.exists(faidx_file):
        cmd = f'samtools faidx {genome}'
        p = Popen(cmd, shell=True)
        return_code = p.wait()
        if return_code != 0:
            print(f'Error: {cmd}')
            sys.exit(1)
    with open(faidx_file, 'r') as f:
        for line in f.readlines():
            lines = line.strip().split('\t')
            result += int(lines[1])
    # 使用 k-mer 工具来模拟 Xbps 长读取到目标基因组的映射，并找到理想的有效基因组大小。然而，通常从整个基因组中去掉简单的重复序列和 Ns，就可以得到一个近似的有效基因组大小。数字上的细微差别不会导致峰值调用的巨大差异，因为这个数字用于估计全基因组噪声水平，与 MACS 建模的局部偏差相比，这个噪声水平通常是最不显著的。
    result = int(result * 0.9)
    return result


def argv_init():
    arg = argparse.ArgumentParser(description='Centromere Find')
    arg.add_argument('-g', '--genome', help='genome path(bwa path)')
    arg.add_argument('-c', '--ckreads1', help='CK reads 1')
    arg.add_argument('-C', '--ckreads2', help='CK reads 2')
    arg.add_argument('-i', '--treatreads1', help='Treat reads 1')
    arg.add_argument('-I', '--treatreads2', help='Treat reads 2')
    arg.add_argument('-t', '--threads', help='threads', default=16)
    arg.add_argument('-f', '--flag', help='flag of genome', default="sfs")
    arg.add_argument(
        '-v', '--view', help='type of view(cent|other)', default="cent")
    args = arg.parse_args()
    return args


def main():
    print("Chip Find Peak START!")
    args = argv_init()
    genome = args.genome
    ckreads1 = args.ckreads1
    ckreads2 = args.ckreads2
    treatreads1 = args.treatreads1
    treatreads2 = args.treatreads2
    threads = int(args.threads)
    flag = args.flag
    view = args.view
    start_time = time.time()
    genome_size = get_genome_size(genome)
    mk_dir(FASTP)
    mk_dir(MACS)
    mk_dir(BAM)
    # 检查reads是否存在
    reads_list = [ckreads1, ckreads2, treatreads1, treatreads2]
    for read in reads_list:
        if not os.path.exists(read):
            print(f'Error: {read} not exists')
            sys.exit(1)
    input_reads1_list = [ckreads1, treatreads1]
    input_reads2_list = [ckreads2, treatreads2]
    output_reads1_list = [os.path.basename(ckreads1).split(
        '_')[0] + '_clean_1.fq.gz', os.path.basename(treatreads1).split('_')[0] + '_clean_1.fq.gz']
    output_reads2_list = [os.path.basename(ckreads2).split(
        '_')[0] + '_clean_2.fq.gz', os.path.basename(treatreads2).split('_')[0] + '_clean_2.fq.gz']
    ck_bam_list = []
    treat_bam_list = []
    ck_bam_list.append(
        flag + '_' + os.path.basename(ckreads1).split('_')[0] + '_rmdup.bam')
    treat_bam_list.append(flag + '_' + os.path.basename(
        treatreads1).split('_')[0] + '_rmdup.bam')
    # fastp
    fastp_theads = 8
    parallel = cacu_threads(threads, fastp_theads)
    with ProcessPoolExecutor(max_workers=parallel) as executor:
        for i in range(len(input_reads1_list)):
            if not check_file(FASTP + "/" + output_reads1_list[i]):
                executor.submit(run_fastp, input_reads1_list[i], input_reads2_list[i], output_reads1_list[i],
                                output_reads2_list[i], input_reads1_list[i].split('_')[0], fastp_theads)
            else:
                print("fastp done:"+input_reads1_list[i])
    print("fastp done")
    output_bam = [flag + '_' +
                  i.split('_')[0] + '_F20.bam' for i in output_reads1_list]
    bwa_theads = 16
    parallel = cacu_threads(threads, bwa_theads)
    with ProcessPoolExecutor(max_workers=parallel) as executor:
        for i in range(len(output_bam)):
            if not check_file(BAM + "/" + output_bam[i]):
                executor.submit(
                    run_bwa, genome, output_reads1_list[i], output_reads2_list[i], output_bam[i].split('_F20')[0], bwa_theads)
            else:
                print("bwa done:"+output_bam[i])
    print("bwa done")
    # dup
    output_rmdup_bam = [flag + '_' +
                        i.split('_')[0] + '_rmdup.bam' for i in output_reads1_list]
    input_F20_bam = [flag + '_' +
                     i.split('_')[0] + '_F20.bam' for i in output_reads1_list]
    dup_bam_threads = 16
    parallel = cacu_threads(threads, dup_bam_threads)
    with ProcessPoolExecutor(max_workers=parallel) as executor:
        for i in range(len(input_F20_bam)):
            if not check_file(BAM + "/" + output_rmdup_bam[i]):
                executor.submit(
                    run_rmdup, input_F20_bam[i], output_rmdup_bam[i])
            else:
                print("rmdup done:"+output_rmdup_bam[i])
    print("rmdup done")
    # macs2
    with ProcessPoolExecutor(max_workers=1) as executor:
        if not check_file(MACS + "/" + flag + "_treat_pileup.bdg"):
            executor.submit(
                run_macs2, ck_bam_list[0], treat_bam_list[0], flag, genome_size)

    ##
    print("macs2 done")
    input_bdg = MACS + "/" + flag + "_chip_sorted.bdg"
    output_wig = input_bdg + ".wig"
    with ProcessPoolExecutor(max_workers=1) as executor:
        if not check_file(output_wig):
            executor.submit(run_bedgraph, input_bdg, genome+".fai", output_wig)
        else:
            print("bedgraph done:"+output_wig)
    print("bedgraph done")

    if view == "cent":
        with ProcessPoolExecutor(max_workers=1) as executor:
            if not check_file(MACS + "/" + flag + "_chip_sorted.bin"):
                executor.submit(chip_plot, flag, genome + ".fai")
            else:
                print("chip convert bin done:"+MACS +
                      "/" + flag + "_chip_sorted.bin")
        print("convert bin done")
    # # 删除BAM文件夹下所有文件
    # # 列出文件
    # file_list = os.listdir(BAM)
    # # 删除
    # for file in file_list:
    #     file_path = os.path.join(BAM, file)
    #     try:
    #         if os.path.isfile(file_path):
    #             os.unlink(file_path)
    #     except Exception as e:
    #         print(e)
    print("time cost:", time.time() - start_time)


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print("User interrupts me! ;-) See you!")
