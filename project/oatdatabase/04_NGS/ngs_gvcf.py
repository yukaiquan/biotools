#!/usr/bin/env python3
import sys
import os
import time
import glob
from subprocess import Popen, PIPE
from concurrent.futures import ProcessPoolExecutor


GATK = '/work/home/acfaa2ssz7/soft/gatk4.py'


def mk_dir(dir) -> str:
    if not os.path.exists(dir):
        os.makedirs(dir)
    return dir


def cacu_threads(threads: int, intheads: int) -> int:
    if threads <= intheads:
        parallel = 1
        intheads = threads
    elif threads > intheads and threads % intheads == 0:
        parallel = threads // intheads
    else:
        parallel = threads // intheads + 1
    return parallel


def run_fastp(fastq1, fastq2, outdir, threads):
    cmd = f'fastp -i {fastq1} -I {fastq2} -o {outdir}/fastp_out1.fq.gz -O {outdir}/fastp_out2.fq.gz -h {outdir}/fastp.html -j {outdir}/fastp.json -w {threads}'
    p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    stdout, stderr = p.communicate()
    print(stdout)
    print(stderr)
    return stdout, stderr


def run_bwa(ref: str, fastq1: str, fastq2: str, output: str, flag: str, threads: int) -> int:
    sample_title = '\'@RG\\tID:' + flag + '\\tSM:' + \
        flag + '\\tLB:' + flag + '\\tPL:Illumina\''
    cmd = f'bwa mem {ref} {fastq1} {fastq2} -t {threads} -R {sample_title} | samtools sort -@ {threads} -o {output}'
    p = Popen(cmd, shell=True)
    return_code = p.wait()
    print(cmd)
    if return_code != 0:
        print(f'Error: {cmd}')
        sys.exit(1)
    return return_code


def rm_dup(input, output, metrics):
    '''
    rm_dup之前的bam需要排过序的
    线程4个就足够了,其实它是单核心的
    '''
    cmd = f'gatk MarkDuplicates -I {input} -O {output} -M {metrics}'
    p = Popen(cmd, shell=True)
    return_code = p.wait()
    return return_code


def gatk4(ref: str, input_bam: str, output: str, threads: int) -> int:
    cmd = f'python {GATK} -P {threads} -I {input_bam} -R {ref} -O {output}'
    p = Popen(cmd, shell=True)
    return_code = p.wait()  # 等待子进程结束，并返回状态码；0表示正常结束，非0表示异常结束
    return return_code


def gatk4_merge(chrom_gvcf, output):
    cmd = f'gatk GatherVcfs -RI TRUE -I {" -I ".join(chrom_gvcf)} -O {output}'
    p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    return_code = p.wait()
    return return_code


def samtools_index(bam: str) -> int:
    '''
    samtools index -c bam.file
    '''
    cmd = f'samtools index -c {bam}'
    p = Popen(cmd, shell=True)
    return_code = p.wait()
    if return_code != 0:
        print(f'Error: {cmd}')
        sys.exit(1)
    return return_code


def g_vcf_index(gvcf) -> int:
    '''
    bcftools version is 1.3.1
    '''
    cmd = f'bcftools index -t {gvcf}'
    if os.path.exists(f'{gvcf}.tbi'):
        return True
    p = Popen(cmd, shell=True)
    return_code = p.wait()
    return return_code


def geno_type_GVCFs(ref, input_gvcf, output_vcf):
    cmd = f'gatk GenotypeGVCFs -OVI False -R {ref} -V {input_gvcf} -O {output_vcf}'
    print(cmd)
    if os.path.exists(output_vcf):
        return True
    p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    stdout, stderr = p.communicate()
    return_code = p.wait()
    if return_code != 0:
        print(f'Error: {cmd} {stderr}')
        sys.exit(1)
    return True


def main():
    ref = '/public/home/acfaa2ssz7/genome/bwaindex/sfs/ACD.sativa.fasta'
    fastq1_list: list = ['/public/home/acfaa2ssz7/YKQNGS/ogle/01_fastq/0002.out.SRR15096497_1.fastq.gz',
                         '/public/home/acfaa2ssz7/YKQNGS/ogle/01_fastq/0003.out.SRR15096497_1.fastq.gz', '/public/home/acfaa2ssz7/YKQNGS/ogle/01_fastq/0004.out.SRR15096497_1.fastq.gz',
                         '/public/home/acfaa2ssz7/YKQNGS/ogle/01_fastq/0005.out.SRR15096497_1.fastq.gz', '/public/home/acfaa2ssz7/YKQNGS/ogle/01_fastq/0006.out.SRR15096497_1.fastq.gz',
                         '/public/home/acfaa2ssz7/YKQNGS/ogle/01_fastq/0007.out.SRR15096497_1.fastq.gz', '/public/home/acfaa2ssz7/YKQNGS/ogle/01_fastq/0008.out.SRR15096497_1.fastq.gz',
                         '/public/home/acfaa2ssz7/YKQNGS/ogle/01_fastq/0009.out.SRR15096497_1.fastq.gz', '/public/home/acfaa2ssz7/YKQNGS/ogle/01_fastq/0010.out.SRR15096497_1.fastq.gz', ]
    fastq2_list: list = ['/public/home/acfaa2ssz7/YKQNGS/ogle/01_fastq/0002.out.SRR15096497_2.fastq.gz',
                         '/public/home/acfaa2ssz7/YKQNGS/ogle/01_fastq/0003.out.SRR15096497_2.fastq.gz', '/public/home/acfaa2ssz7/YKQNGS/ogle/01_fastq/0004.out.SRR15096497_2.fastq.gz',
                         '/public/home/acfaa2ssz7/YKQNGS/ogle/01_fastq/0005.out.SRR15096497_2.fastq.gz', '/public/home/acfaa2ssz7/YKQNGS/ogle/01_fastq/0006.out.SRR15096497_2.fastq.gz',
                         '/public/home/acfaa2ssz7/YKQNGS/ogle/01_fastq/0007.out.SRR15096497_2.fastq.gz', '/public/home/acfaa2ssz7/YKQNGS/ogle/01_fastq/0008.out.SRR15096497_2.fastq.gz',
                         '/public/home/acfaa2ssz7/YKQNGS/ogle/01_fastq/0009.out.SRR15096497_2.fastq.gz', '/public/home/acfaa2ssz7/YKQNGS/ogle/01_fastq/0010.out.SRR15096497_2.fastq.gz', ]
    sample_name = 'ogle'
    bam_out = '02_bam'
    dup_out = '03_dup'
    g_vcf = '04_gvcf'
    threads = 32
    bam_list: list = []
    bam_name: list = []
    dup_bam_list: list = []
    gvcf_list: list = []
    # 创建输出目录
    print('Create output directory', mk_dir(bam_out))
    bwa_theads = 16
    parallel = cacu_threads(threads, bwa_theads)
    with ProcessPoolExecutor(max_workers=parallel) as executor:
        for i in range(len(fastq1_list)):
            fastq1 = fastq1_list[i]
            fastq2 = fastq2_list[i]
            # 提取fastq1 和fastq2 共有部分作为文件名
            sample1_name = os.path.basename(fastq1).split('_')[0]
            sample2_name = os.path.basename(fastq2).split('_')[0]
            # 判断两个文件名是否一致
            if sample1_name != sample2_name:
                print('Error: fastq1 and fastq2 name is not equal, please check !')
                sys.exit(1)
            output_file = os.path.join(bam_out, sample1_name + '.bam')
            bam_list.append(output_file)
            bam_name.append(sample1_name + '.bam')
            executor.submit(run_bwa, ref, fastq1, fastq2,
                            sample1_name + '.bam', sample_name, bwa_theads)
    print('#'*25, 'bwa done !', '#'*25)
    # 移动文件到bam_out目录
    for bam in bam_name:
        os.system('mv %s %s' % (bam, bam_out))
    # 创建index
    index_threads = 1
    parallel = cacu_threads(threads, index_threads)
    with ProcessPoolExecutor(max_workers=parallel) as executor:
        for bam in bam_list:
            executor.submit(samtools_index, bam)
    print('#'*25, 'samtools index done !', '#'*25)
    # bam to dup
    print('Create output directory', mk_dir(dup_out))
    dup_theads = 4
    parallel = cacu_threads(threads, dup_theads)
    with ProcessPoolExecutor(max_workers=parallel) as executor:
        for i in range(len(bam_list)):
            bam = bam_list[i]
            output_file = os.path.join(
                dup_out, os.path.basename(bam) + '.dup.bam')
            dup_bam_list.append(output_file)
            executor.submit(rm_dup, bam, output_file, dup_theads)
    print('#'*25, 'rm_dup done !', '#'*25)
    # 创建index
    index_threads = 1
    parallel = cacu_threads(threads, index_threads)
    with ProcessPoolExecutor(max_workers=parallel) as executor:
        for bam in dup_bam_list:
            executor.submit(samtools_index, bam)
    # dup to gvcf
    print('Create output directory', mk_dir(g_vcf))
    g_vcf_threads = 16
    parallel = cacu_threads(threads, g_vcf_threads)
    with ProcessPoolExecutor(max_workers=parallel) as executor:
        for i in range(len(dup_bam_list)):
            bam = dup_bam_list[i]
            output_file = os.path.basename(bam)
            gvcf_list.append(output_file)
            executor.submit(gatk4, ref, bam, output_file, g_vcf_threads)
    # 检测生成的gvcf文件大小
    for gvcf in gvcf_list:
        chr_gvcf_list = glob.glob(os.path.join(g_vcf, gvcf + '*.g.vcf.gz'))
        for chr_gvcf in chr_gvcf_list:
            # 小于1M
            if os.path.getsize(chr_gvcf) < 1024 * 1024:
                print('Error: gvcf file size is less than 1M, please check !', chr_gvcf)
            # 删除染色体小文件
            # os.remove(chr_gvcf)
            print('gvcf file size is ok !', chr_gvcf)
    # 所有任务结束
    print('#'*25, 'all done !', '#'*25)


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print('KeyboardInterrupt')
        sys.exit(0)
