#!/usr/bin/env python3
import sys
import os
from subprocess import Popen, PIPE
from concurrent.futures import ProcessPoolExecutor
import argparse


# samtools 版本: 1.12 太低了不支持bwa mem的输出顺便排序

GATK: str = '/public/home/acfaa2ssz7/soft/gatk4.py'
CHR_LIST: list = ['chr1A', 'chr2A', 'chr3A', 'chr4A', 'chr5A', 'chr6A', 'chr7A',
                  'chr1C', 'chr2C', 'chr3C', 'chr4C', 'chr5C', 'chr6C', 'chr7C',
                  'chr1D', 'chr2D', 'chr3D', 'chr4D', 'chr5D', 'chr6D', 'chr7D', 'chrUn']


def mk_dir(dir) -> str:
    if not os.path.exists(dir):
        os.makedirs(dir)
    return dir


def mv_file(file, outdir) -> str:
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    os.system(f'mv {file} {outdir}')
    return outdir


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


def run_fastp(fastq1, fastq2, outdir, out_name1, out_name2, split, threads) -> bool:
    if split == 1:
        cmd = f'fastp -i {fastq1} -I {fastq2} -o {outdir}/{out_name1} -O {outdir}/{out_name2} -h {outdir}/{out_name1}_{out_name2}_fastp.html -j {outdir}/{out_name1}_{out_name2}_fastp.json -w {threads}'
    else:
        cmd = f'fastp -i {fastq1} -I {fastq2} -o {outdir}/{out_name1} -O {outdir}/{out_name2} -h {outdir}/{out_name1}_{out_name2}_fastp.html -j {outdir}/{out_name1}_{out_name2}_fastp.json -s {split} -w {threads}'
    p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    stdout, stderr = p.communicate()
    return_code = p.wait()
    print(stdout)
    print(stderr)
    if return_code != 0:
        print(f'Error: {cmd}')
        sys.exit(1)
    else:
        print(f'Fastp finished: {cmd}')
        return True


def run_bwa(ref: str, fastq1: str, fastq2: str, output: str, flag: str, threads: int) -> int:
    sample_title = '\'@RG\\tID:' + flag + '\\tSM:' + \
        flag + '\\tLB:' + flag + '\\tPL:Illumina\''
    cmd = f'bwa mem {ref} {fastq1} {fastq2} -t {threads} -R {sample_title} | samtools sort -@ {threads} -o {output}'
    p = Popen(cmd, shell=True)
    return_code = p.wait()
    stdout, stderr = p.communicate()
    print(cmd)
    if return_code != 0:
        print(f'Error: {cmd}')
        print(stderr)
        sys.exit(1)
    return True


def rm_dup(input, output, metrics):
    '''
    rm_dup之前的bam需要排过序的
    线程8个就足够了,其实它是单核心的 但是会占用很多内存 有时候CPU又会被占用
    '''
    cmd = f'gatk --java-options -Xmx12G MarkDuplicates -I {input} -O {output} -M {metrics}'
    p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    stdout, stderr = p.communicate()
    return_code = p.wait()
    if return_code != 0:
        print(f'Error: {cmd} {stdout}{stderr}')
        sys.exit(1)
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


def samtools_index(bam: str) -> True:
    '''
    samtools index -c bam.file
    '''
    if os.path.exists(f'{bam}.csi'):
        return True
    cmd = f'samtools index -c {bam}'
    p = Popen(cmd, shell=True)
    return_code = p.wait()
    if return_code != 0:
        print(f'Error: {cmd}')
        sys.exit(1)
    return True


def g_vcf_index(gvcf: str) -> int:
    '''
    bcftools version is 1.3.1
    '''
    cmd = f'bcftools index -t {gvcf}'
    if os.path.exists(f'{gvcf}.tbi'):
        return True
    p = Popen(cmd, shell=True)
    return_code = p.wait()
    return return_code


def geno_type_GVCFs(ref: str, input_gvcf: str, output_vcf: str) -> bool:
    cmd = f'gatk GenotypeGVCFs -OVI False -R {ref} -V {input_gvcf} -O {output_vcf}'
    print(cmd)
    if os.path.exists(output_vcf):
        return True
    p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    stdout, stderr = p.communicate()
    print(stdout)
    return_code = p.wait()
    if return_code != 0:
        print(f'Error: {cmd} {stderr}')
        sys.exit(1)
    return True

# 写一个函数将多个染色体的运行命令打包为list再利用多进程的方式将所有命令提交进去，充分利用计算资源


def genetion_gvcf(ref: str, input_bam: str, output_gvcf: str, chrom: str) -> str:
    out_put_name = output_gvcf + '.' + chrom + '.g.vcf.gz'
    cmd = f'gatk --java-options -Xmx6G HaplotypeCaller -I {input_bam} -O {out_put_name} -R {ref} --emit-ref-confidence GVCF -OVI False -L {chrom}'
    # print(cmd)
    return cmd


def gvcf_cmd_list(ref, input_bam, output_gvcf) -> list:
    cmd_list = []
    for chrom in CHR_LIST:
        cmd = genetion_gvcf(ref, input_bam, output_gvcf, chrom)
        cmd_list.append(cmd)
    return cmd_list


def pysh(cmd) -> bool:
    '''
    执行输入的shell命令,并返回执行状态bool
    '''
    p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    stdout, stderr = p.communicate()
    return_code = p.wait()
    if return_code != 0:
        print(f'Error: {cmd} {stderr}')
        sys.exit(1)
    return True


def gatk4_merge(chrom_gvcf: list, output: str) -> int:
    cmd = f'gatk GatherVcfs -RI TRUE -I {" -I ".join(chrom_gvcf)} -O {output}'
    p = Popen(cmd, shell=True)
    return_code = p.wait()
    return return_code


def sample_fasq_list(fastq1_list: list, fastq2_list: list, flag_list: list, sample_input_fq: str) -> list:
    with open(sample_input_fq, 'r') as f:
        for line in f:
            fastq1_list.append(line.strip().split()[0])
            fastq2_list.append(line.strip().split()[1])
            flag_list.append(line.strip().split()[2])
    return fastq1_list, fastq2_list, flag_list


def file_exists(file: str) -> bool:
    if os.path.exists(file):
        return True
    else:
        return False


def main():
    # 读取命令行参数
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='input fastq file list')
    parser.add_argument('-r', '--ref', help='reference genome')
    parser.add_argument('-t', '--threads', help='threads', default=32)
    parser.add_argument(
        '-s', '--split', help='split fastq default 1 is not split fastq', default=1)
    args = parser.parse_args()
    # ref = '/mnt/disk3/new_genome/ACD/bwaindex/ACD.sativa.fasta'
    # ref = '/public/home/acfaa2ssz7/genome/bwaindex/sfs/ACD.sativa.fasta'
    ref = args.ref
    # sample_input_fq = './sample_input_fq.txt'
    sample_input_fq = args.input
    fastq1_list = []
    fastq2_list = []
    flag_list = []
    out_fastp1_list = []
    out_fastp2_list = []
    fastq1_list, fastq2_list, flag_list = sample_fasq_list(
        fastq1_list, fastq2_list, flag_list, sample_input_fq)
    # fastq1_list: list = [
    #     '/public/home/acfaa2ssz7/YKQNGS/ogle/SRR15096498/01_fastq/0004.out.SRR15096498_1.fastq.gz',
    #     '/public/home/acfaa2ssz7/YKQNGS/ogle/SRR15096498/01_fastq/0009.out.SRR15096498_1.fastq.gz']
    # fastq2_list: list = ['/public/home/acfaa2ssz7/YKQNGS/ogle/SRR15096498/01_fastq/0004.out.SRR15096498_2.fastq.gz',
    #                      '/public/home/acfaa2ssz7/YKQNGS/ogle/SRR15096498/01_fastq/0009.out.SRR15096498_2.fastq.gz', ]
    # 检测输入文件是否存在
    if file_exists(ref) is False:
        print('Reference file does not exist')
        sys.exit(1)
    # sample_name = 'sang'
    fastq_out = '01_fastq'
    bam_out = '02_bam'
    dup_out = '03_dup'
    g_vcf = '04_gvcf'
    # threads = 32
    threads = int(args.threads)
    # 决定用不用fastp拆分输入文件 1不拆 多个为拆分的数量
    # is_split = 1
    is_split = int(args.split)
    bam_list: list = []
    bam_name: list = []
    dup_bam_list: list = []
    gvcf_list: list = []
    # 创建输出目录
    fastp_thread = 16
    print('Create output directory', mk_dir(fastq_out))
    parallel = cacu_threads(threads, fastp_thread)
    with ProcessPoolExecutor(max_workers=parallel) as executor:
        for i in range(len(fastq1_list)):
            fastq1 = fastq1_list[i]
            fastq2 = fastq2_list[i]
            # 提取fastq1 和fastq2 共有部分作为文件名
            sample1_name = os.path.basename(fastq1).split('_')[0]
            sample2_name = os.path.basename(fastq2).split('_')[0]
            # 判断两个文件名是否一致
            if sample1_name != sample2_name:
                print('Error: fastq1 and fastq2 name is not equal')
                sys.exit(1)
            # 创建输出文件名
            fastq1_out = sample1_name + '_out_1.fastq.gz'
            fastq2_out = sample2_name + '_out_2.fastq.gz'
            # 将输出文件名添加到列表中
            out_fastp1_list.append(os.path.join(fastq_out, fastq1_out))
            out_fastp2_list.append(os.path.join(fastq_out, fastq2_out))
            # # 创建fastp命令
            if executor.submit(run_fastp, fastq1, fastq2, fastq_out,
                               fastq1_out, fastq2_out, is_split, fastp_thread):
                print('fastp is running')
            else:
                print('fastp is error')
                sys.exit(1)
    parallel = cacu_threads(threads, fastp_thread)
    # 创建输出目录
    print('Create output directory', mk_dir(bam_out))
    bwa_theads = 32
    parallel = cacu_threads(threads, bwa_theads)
    with ProcessPoolExecutor(max_workers=parallel) as executor:
        for i in range(len(out_fastp1_list)):
            fastq1 = out_fastp1_list[i]
            fastq2 = out_fastp2_list[i]
            sample_name = flag_list[i]
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
            if executor.submit(run_bwa, ref, fastq1, fastq2,
                               output_file, sample_name, bwa_theads):
                print('Bwa mapping done !')
            else:
                print('Error: Bwa mapping failed !')
                sys.exit(1)
    print('#'*25, 'bwa done !', '#'*25)
    # # 移动文件到bam_out目录
    # for bam in bam_name:
    #     os.system('mv %s %s' % (bam, bam_out))
    #     print('move bam file: ', bam, 'done !')
    # 创建index 线程数为1
    index_threads = 1
    parallel = cacu_threads(threads, index_threads)
    with ProcessPoolExecutor(max_workers=parallel) as executor:
        for bam in bam_list:
            # executor.submit(samtools_index, bam)
            print('index bam file: ', bam, 'done !')
    print('#'*25, 'samtools index done !', '#'*25)
    # bam to dup threads 2
    print('Create output directory', mk_dir(dup_out))
    dup_theads = 16
    parallel = cacu_threads(threads, dup_theads)
    with ProcessPoolExecutor(max_workers=parallel) as executor:
        for bam in bam_list:
            output_file = os.path.join(
                dup_out, os.path.basename(bam) + '.dup.bam')
            output_matrix = os.path.join(
                dup_out, os.path.basename(bam) + '.dup.matrix')
            dup_bam_list.append(output_file)
            executor.submit(rm_dup, bam, output_file, output_matrix)
    print('#'*25, 'rm_dup done !', '#'*25)
    # 创建index threads 1
    index_threads = 1
    parallel = cacu_threads(threads, index_threads)
    with ProcessPoolExecutor(max_workers=parallel) as executor:
        for bam in dup_bam_list:
            executor.submit(samtools_index, bam)
            print('index bam file: ', bam, 'done !')
    # dup to gvcf threads 1
    print('Create output directory', mk_dir(g_vcf))
    g_vcf_threads = 1
    # dup_bam_list = [
    #     '/public/home/acfaa2ssz7/YKQNGS/ogle/03_dup/0010.out.SRR15096498.bam.dup.bam']
    # 存储对应样本生成染色体gvcf的命令
    gvcf_dict: dict = {}
    # 存储对应样本对应生成的gvcf文件
    gvcf_chr: dict = {}
    for bam in dup_bam_list:
        output_name = os.path.basename(bam).split('.')[0]
        gvcf_dict[output_name] = gvcf_cmd_list(ref, bam, output_name)
        output_name_list = [output_name + '.' + chr +
                            '.g.vcf.gz' for chr in CHR_LIST]
        gvcf_chr[output_name] = output_name_list
    parallel = cacu_threads(threads, g_vcf_threads)
    with ProcessPoolExecutor(max_workers=parallel) as executor:
        # 迭代字典里所有值
        for gvcf in gvcf_dict.values():
            for cmd in gvcf:
                executor.submit(pysh, cmd)
                print('#'*25, cmd, 'gvcf done !', '#'*25)
    gvcf_merge_threads = 2
    parallel = cacu_threads(threads, gvcf_merge_threads)
    with ProcessPoolExecutor(max_workers=parallel) as executor:
        # 按照字典的键值对迭代
        for gvcf in gvcf_chr.keys():
            print('#'*25, gvcf, 'merge start !', '#'*25)
            # print(gvcf_chr[gvcf])
            executor.submit(gatk4_merge, gvcf_chr[gvcf], gvcf + '.g.vcf.gz')
            print('#'*25, gvcf, 'merge done !', '#'*25)
    # 移动文件到gvcf_out目录
    for gvcf in gvcf_chr.keys():
        # 判断文件是否存在
        if os.path.exists(gvcf + '.g.vcf.gz'):
            mv_file(gvcf + '.g.vcf.gz', g_vcf)
            print('#'*25, gvcf, 'move done !', '#'*25)
        # 判断值的文件是否存在
        for chr_gvcf in gvcf_chr[gvcf]:
            if os.path.exists(chr_gvcf):
                mv_file(chr_gvcf, g_vcf)
                print('#'*25, chr_gvcf, 'move done !', '#'*25)
    # 检测生成的gvcf文件大小
    for gvcf in gvcf_chr.keys():
        n: int = 0
        for chr_single in gvcf_chr[gvcf]:
            chr_gvcf = os.path.join(g_vcf, chr_single)
            if os.path.exists(chr_gvcf):
                if os.path.getsize(chr_gvcf) < 1024 * 10:
                    print(
                        'Error: gvcf file size is less than 10k, please check !', chr_gvcf)
                else:
                    # 删除染色体小文件
                    os.remove(chr_gvcf)
                    # print(os.path.getsize(chr_gvcf))
                    print('gvcf file size is ok !', chr_gvcf)
                    n += 1
        if n == len(gvcf_chr[gvcf]):
            mv_file(gvcf + '.g.vcf.gz', g_vcf)
            print('#'*25, gvcf, 'All size done !', '#'*25)
            # 删除染色体小文件
            for chr_single in gvcf_chr[gvcf]:
                chr_gvcf = os.path.join(g_vcf, chr_single)
                # 判断文件是否存在
                if os.path.exists(chr_gvcf):
                    # 删除染色体小文件
                    os.remove(chr_gvcf)
                print('remove file: ', chr_gvcf)
        else:
            print('Error: gvcf file size is less than 1M, please check !', gvcf)
    # 创建索引
    tbi_threads = 1
    parallel = cacu_threads(threads, tbi_threads)
    with ProcessPoolExecutor(max_workers=parallel) as executor:
        # 按照字典的键值对迭代
        for gvcf in gvcf_chr.keys():
            # 判断文件是否存在
            if os.path.exists(os.path.join(g_vcf, gvcf + '.g.vcf.gz.tbi')):
                print('#'*25, gvcf, 'index end !', '#'*25)
            else:
                print('Error: gvcf tbi file exists !', gvcf)
                executor.submit(g_vcf_index, os.path.join(
                    g_vcf, gvcf + '.g.vcf.gz'))
    # 合并所有gvcf文件
    # gvcf_list = [os.path.join(g_vcf, gvcf + '.g.vcf.gz') for gvcf in gvcf_chr.keys()]
    # # dele dup bam文件
    for bam in dup_bam_list:
        # 判断文件是否存在
        if os.path.exists(bam):
            # 删除dup bam文件
            os.remove(bam)
            print('remove file: ', bam)
    # 所有任务结束
    print('#'*25, 'all done !', '#'*25)


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print('KeyboardInterrupt')
        sys.exit(0)
