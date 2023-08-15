#!/usr/bin/env python3

import sys
import os
from subprocess import Popen, PIPE
# 多线程
from concurrent.futures import ProcessPoolExecutor
import shutil
import time

FASTP_OUT_DIR = '01_fastp'
BAM_OUT_DIR = '02_bam'
FEATURE_COUNTS = '03_feature_counts'
MATRIX = '04_matrix'
DEG = '05_deg'
DUP_BAM = '06_dup_bam'
ADD_HEADER = '07_add_header'
GVCF = '08_gvcf'
VCF = '09_vcf'

R_base = '/public/workspace/kqyu/anaconda3/envs/yukaiquan/bin/Rscript'
Feature_counts = '/mnt/disk3/soft/RNAseqScript/run-featurecounts.R'
Merge_featurecounts = '/mnt/disk3/soft/RNAseqScript/abundance_estimates_to_matrix.pl'
DEG = '/mnt/disk3/soft/RNAseqScript/run_DE_analysis.pl'
GATK = '/mnt/disk3/soft/RNAseqScript/gatk4.py'
Tassle ='/mnt/disk3/soft/tassel-5.2.40-3/run_pipeline.pl'
BCFTOOLS_131 ='/usr/local/bin/bcftools'

# 判断上方文件是否存在
if not os.path.exists(R_base):
    print('Error: R_base not exist')
    sys.exit(1)
if not os.path.exists(Feature_counts):
    print('Error: Feature_counts not exist')
    sys.exit(1)
if not os.path.exists(Merge_featurecounts):
    print('Error: Merge_featurecounts not exist')
    sys.exit(1)
if not os.path.exists(GATK):
    print('Error: GATK not exist')
    sys.exit(1)

import argparse
import subprocess
from concurrent.futures import ProcessPoolExecutor

def get_chr_name(bam:str) -> list:
    """
    # samtools view -H 02_bam/RBSR.sam | grep '^@SQ' | sed -r 's#^@SQ\tSN:(.*)\tLN.*$#\1#g' 获取对应物种的染色体名称
    获取bam文件中的染色体名称
    :param bam: bam文件
    :return: 染色体名称列表
    """
    cmd = 'samtools view -H ' + bam + ' | grep \'^@SQ\' | sed -r \'s#^@SQ\\tSN:(.*)\\tLN.*$#\\1#g\''
    print(cmd)
    chr_name = subprocess.getoutput(cmd).split('\n')
    print(chr_name)
    return chr_name

def gatk4(chrom:str,input:str,output:str,ref:str):
    subprocess.run(['gatk','--java-options','-Xmx8G','HaplotypeCaller','-I',input,'-O',output + '.' + chrom + '.g.vcf.gz','-R',ref,'--emit-ref-confidence','GVCF','-OVI','False', '-L',chrom],shell=False)

def gatk_split(bam:str,ref:str,process:int,output:str):
    """
    :param bam: bam文件
    :param ref: 参考基因组
    :param process: 并行进程数
    :param output: 输出文件名
    """
    chrom_list = get_chr_name(bam)
    # 燕麦的21条染色体 按需求替换名字要对应参考基因组
    chrom_gvcf = []
    for chrom in chrom_list:
        chrom_gvcf.append(output + '.' + chrom + '.g.vcf.gz')
    try:
        with ProcessPoolExecutor(max_workers=process) as pool:
            future1 = pool.map(gatk4,chrom_list,[bam] * len(chrom_list),[output] * len(chrom_list),[ref] * len(chrom_list))
            #print(list(future1))
    except Exception as e: 
        print(e)
    subprocess.run('gatk GatherVcfs -RI TRUE -I ' + ' -I '.join(chrom_gvcf) + ' -O ' + output + '.g.vcf.gz',shell=True)
    #subprocess.run('rm -fr ' + ' '.join(chrom_gvcf),shell=True)

def cacu_threads(threads: int, intheads: int) -> int:
    if threads <= intheads:
        parallel = 1
        intheads = threads
    elif threads > intheads and threads % intheads == 0:
        parallel = threads // intheads
    else:
        parallel = threads // intheads + 1
    return parallel


def fastp(fastq_1: str, fastq_2: str, threads: int) -> list:
    cmd = 'fastp -i ' + fastq_1 + ' -I ' + fastq_2 + ' -o ' + FASTP_OUT_DIR + '/' + fastq_1.split('/')[-1].split('_')[0] + 'clean_1.fastq.gz -O ' + FASTP_OUT_DIR + '/' + fastq_2.split('_')[-1].split(
        '.')[0] + 'clean_2.fastq.gz -h ' + FASTP_OUT_DIR + '/' + fastq_1.split('/')[-1].split('_')[0] + '.html -j ' + FASTP_OUT_DIR + '/' + fastq_1.split('/')[-1].split('_')[0] + '.json' + ' -w ' + str(threads)
    print(cmd)
    output = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    stdout, stderr = output.communicate()
    if output.returncode != 0:
        print(stderr)
        sys.exit(1)
    output.wait()
    return [FASTP_OUT_DIR + '/' + fastq_1.split('/')[-1].split('.')[0] + '_clean.fastq.gz', FASTP_OUT_DIR + '/' + fastq_2.split('/')[-1].split('.')[0] + '_clean.fastq.gz']


def hisat2(fastq_1: str, fastq_2: str, threads: int, genome: str, output: str) -> str:
    cmd = 'hisat2 -p ' + str(threads) + ' --dta -x ' + genome + \
        ' -1 ' + fastq_1 + ' -2 ' + fastq_2 + ' -S ' + output
    print(cmd)
    out = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    stdout, stderr = out.communicate()
    if out.returncode != 0:
        print(stderr)
        sys.exit(1)
    out.wait()
    return output


def samtools_sort(bam: str, threads: int) -> str:
    cmd = 'samtools sort -@ ' + \
        str(threads) + ' -o ' + bam.split('.')[0] + '_sorted.bam ' + bam
    print(cmd)
    out = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    stdout, stderr = out.communicate()
    if out.returncode != 0:
        print(stderr)
        sys.exit(1)
    out.wait()
    return bam.split('.')[0] + '_sorted.bam'


def gatk_remove_duplicates(bam, out):
    cmd = f'gatk MarkDuplicates -I {bam} -O {out} -M {out}.metrics.txt'
    p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    stdout, stderr = p.communicate()
    p.wait()  # 等待子进程结束，并返回状态码；
    return stdout, stderr


def samtools_index(bam: str) -> str:
    cmd = 'samtools index -c ' + bam
    print(cmd)
    out = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    stdout, stderr = out.communicate()
    if out.returncode != 0:
        print(stderr)
        sys.exit(1)
    out.wait()
    return bam + '.csi'

def gatk_add_read_groups(bam:str,out:str,header:str):
    cmd = f'gatk AddOrReplaceReadGroups -I {bam} -O {out} -SO coordinate -ID {header} -LB {header} -PL illumina -PU {header} -SM {header}'
    p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    stdout, stderr = p.communicate()
    p.wait()  # 等待子进程结束，并返回状态码；
    return stdout, stderr

def gatk_genotype_vcf(gvcf:str,out:str,ref:str):
    cmd = f'gatk GenotypeGVCFs -OVI False -R {ref} -V {gvcf} -O {out}'
    print(cmd)
    p = Popen(cmd, shell=True)
    p.wait()  # 等待子进程结束

def mkdir(dir: str):
    if not os.path.exists(dir):
        os.mkdir(dir)


def mk_samples_dict(sequence_list: str) -> dict:
    samples_dict = {}
    with open(sequence_list, 'r') as f:
        for line in f:
            line = line.strip().split()
            samples_dict[line[0]] = line[1:]
    return samples_dict

def run_featurecounts(R_base, bam, gtf, name, threads):
    cmd = f'{R_base} {Feature_counts} -b {bam} -g {gtf} -o {name} -t {threads}'
    p = Popen(cmd, shell=True)
    return_code = p.wait()  # 等待子进程结束，并返回状态码；
    return return_code

def run_merge_featurecounts(out, sample_list) -> list:
    # perl ~/soft/abundance_estimates_to_matrix.pl --est_method featureCounts --cross_sample_norm TMM --out_prefix dw6.matrix --quant_files
    cmd = f'perl {Merge_featurecounts} --est_method featureCounts --cross_sample_norm TMM --out_prefix {out} --quant_files {sample_list}'
    p = Popen(cmd, shell=True)
    return_code = p.wait()  # 等待子进程结束，并返回状态码；
    return [out + '.counts.matrix', out + '.TMM.EXPR.matrix', out + '.TPM.not_cross_norm', out + '.TPM.not_cross_norm.TMM_info.txt']


def list_2_file(list: list, file: str):
    with open(file, 'w') as f:
        for i in list:
            f.write(FEATURE_COUNTS + '/' + i + '\n')

def mv_file(file: str, dir: str):
    cmd = f'mv {file} {dir}'
    p = Popen(cmd, shell=True)
    return_code = p.wait()  # 等待子进程结束，并返回状态码；
    return return_code
def bcfTools_index(bcf: str):
    cmd = BCFTOOLS_131 + ' index -t ' + bcf
    print(cmd)
    out = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    stdout, stderr = out.communicate()
    if out.returncode != 0:
        print(stderr)
        sys.exit(1)
    out.wait()

def gatk_merge_vcf(vcf_list: list,genome:str, out: str):
    cmd = 'gatk CombineGVCFs -OVI False -V ' + ' -V '.join(vcf_list) + ' -R ' + genome + ' -O ' + out 
    print(cmd)
    out = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    stdout, stderr = out.communicate()
    if out.returncode != 0:
        print(stderr)
        sys.exit(1)
    out.wait()

def bcfTools_call(bam: str, out: str, genome: str):
    ## 参数
    # -f: 指定参考基因组
    # -b: bam list的文件，样本较多时可以使用
    # -C: --adjust-MQ 矫正的MQ值，推介50
    # -q: --mim-MQ  MQ质量值
    # -Q: --min-BQ  base质量值
    # -r: --regions call 特定染色体或者区域的变异；如chr1; chr1:100-20000
    # -R: --regionbs-file 当有多个区域，写如一个文件；tab分割： chi  start. end
    # -a: --annoteta 在INFO中添加DP4，AD等信息
    # -O: --output-type  有b, u,z,v
    # ## call 参数
    # -s, --samples : 列出需要call的样本，默认全部
    # -S, --samples-file 支队文件中列出的样本进行call
    # -c : 对应consensus-caller算法
    # -m 对应multiallelic-caller算法，后者更适合多种allel和罕见变异的calling。
    # -v: 只输出变异位点，不输出参考位点
    cmd = 'bcftools mpileup -f ' + genome + ' -Ou ' + bam + ' | bcftools call -mv -Oz -o ' + out
    print(cmd)
    out = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    stdout, stderr = out.communicate()
    if out.returncode != 0:
        print(stderr)
        sys.exit(1)
    out.wait()

def bcfTools_merge(bcf_list: list, out: str):
    cmd = 'bcftools merge -Oz -o ' + out + ' ' + ' '.join(bcf_list)
    print(cmd)
    out = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    stdout, stderr = out.communicate()
    if out.returncode != 0:
        print(stderr)
        sys.exit(1)
    out.wait()

def rm_files(files: list):
    if len(files) == 0:
        return
    for file in files:
        if os.path.exists(file):
            os.remove(file)

def vcf_2_hapmap(vcf: str, out: str):
    # /mnt/disk3/soft/tassel-5.2.40-3/run_pipeline.pl -Xms64G -Xmx64G -fork1 -vcf 09_vcf/all.vcf.gz -export 09_vcf/all  -exportType Hapmap
    cmd = Tassle + ' -Xms64G -Xmx64G -fork1 -vcf ' + vcf + ' -export ' + out + ' -exportType Hapmap'
    print(cmd)
    out = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    stdout, stderr = out.communicate()
    if out.returncode != 0:
        print(stderr)
        sys.exit(1)
    out.wait()


def main():
    print("### bsrSeq.py ###")
    print("### This script will run the BSR algorithm on a single sequence ###")
    print("### Usage: python3 bsrSeq.py <genome> <sequence list> ###")
    print("sequence list format: <uniq name> <fastq_1> <fastq_2> condition1 condition2...")
    time_start = time.time()
    genome = sys.argv[1]
    genome_index = sys.argv[2]
    genome_gtf = sys.argv[3]
    sequence_list = sys.argv[4]
    threads = int(sys.argv[5])
    mkdir(FASTP_OUT_DIR)
    samples_dict = mk_samples_dict(sequence_list)
    #### fastp ####
    print("### fastp ###")
    fastp_dict = {}
    fastp_threads = 16
    parallel = cacu_threads(threads, fastp_threads)
    with ProcessPoolExecutor(max_workers=parallel) as executor:
        for sample in samples_dict:
            fastp_out_1 = FASTP_OUT_DIR + '/' + \
                samples_dict[sample][0].split(
                    '/')[-1].split('_')[0] + 'clean_1.fastq.gz'
            fastp_out_2 = FASTP_OUT_DIR + '/' + \
                samples_dict[sample][1].split(
                    '/')[-1].split('_')[0] + 'clean_2.fastq.gz'
            fastp_dict[sample] = [fastp_out_1, fastp_out_2]
            # 如果结果文件已经存在，跳过
            # print(samples_dict[sample][0].split('/')[-1].split('_')[0]+ '.html')
            if os.path.exists(FASTP_OUT_DIR + '/' +samples_dict[sample][0].split('/')[-1].split('_')[0]+ '.html'):
                print(fastp_out_1 + ' and ' +
                      fastp_out_2 + ' already exists, skip!')
                # 跳过本轮循环
                continue
            else:
                print("### fastp ###")
                executor.submit(
                    fastp, samples_dict[sample][0], samples_dict[sample][1], threads=fastp_threads)
    #### hisat2 ####
    print("### hisat2 ###")
    mkdir(BAM_OUT_DIR)
    hisat2_dict = {}
    hisat2_threads = 16
    parallel = cacu_threads(threads, hisat2_threads)
    with ProcessPoolExecutor(max_workers=parallel) as executor:
        for sample in fastp_dict:
            hisat2_out = BAM_OUT_DIR + '/' + sample + '.sam'
            hisat2_dict[sample] = hisat2_out
            # 如果结果文件已经存在，跳过
            if os.path.exists(BAM_OUT_DIR + '/' + sample + '_sorted.bam'):
                print(hisat2_out + ' already exists, skip!')
                continue
            else:
                executor.submit(hisat2, fastp_dict[sample][0], fastp_dict[sample][1],
                                threads=hisat2_threads, genome=genome_index, output=hisat2_out)
    #### samtools sort ####
    print("### samtools sort ###")
    samtools_sort_dict = {}
    samtools_sort_threads = 8
    parallel = cacu_threads(threads, samtools_sort_threads)
    with ProcessPoolExecutor(max_workers=parallel) as executor:
        for sample in hisat2_dict:
            samtools_sort_out = BAM_OUT_DIR + '/' + sample + '_sorted.bam'
            samtools_sort_dict[sample] = samtools_sort_out
            # 如果结果文件已经存在，跳过
            if os.path.exists(samtools_sort_out):
                print(samtools_sort_out + ' already exists, skip!')
                continue
            else:
                executor.submit(
                    samtools_sort, hisat2_dict[sample], threads=samtools_sort_threads)
    #### samtools index ####
    print("### samtools index ###")
    samtools_index_threads = 1
    parallel = cacu_threads(threads, samtools_index_threads)
    with ProcessPoolExecutor(max_workers=parallel) as executor:
        for sample in samtools_sort_dict:
            samtools_index_out = samtools_sort_dict[sample] + '.csi'
            # 如果结果文件已经存在，跳过
            if os.path.exists(samtools_index_out):
                print(samtools_index_out + ' already exists, skip!')
                continue
            else:
                executor.submit(samtools_index, samtools_sort_dict[sample])
    #### run-featureCounts ####
    print("### run-featureCounts ###")
    mkdir(FEATURE_COUNTS)
    feature_counts_list = []
    featureCounts_threads = 8
    parallel = cacu_threads(threads, featureCounts_threads)
    with ProcessPoolExecutor(max_workers=parallel) as executor:
        for sample in samtools_sort_dict:
            featureCounts_out =  sample
            feature_counts_list.append(featureCounts_out + '.count')
            # 如果结果文件已经存在，跳过
            if os.path.exists(FEATURE_COUNTS + '/' + featureCounts_out + '.count'):
                print(featureCounts_out + ' already exists, skip!')
                continue
            else:
                executor.submit(
                    run_featurecounts, R_base, samtools_sort_dict[sample], genome_gtf, featureCounts_out, threads=featureCounts_threads)
    for file in feature_counts_list:
        if os.path.exists(FEATURE_COUNTS + '/' + file):
            print(file + ' already exists, skip!')
            continue
        else:
            shutil.move(file, FEATURE_COUNTS)
    #### matrix ####
    print("### matrix ###")
    mkdir(MATRIX)
    if os.path.exists(MATRIX + '/RNA.counts.matrix'):
        print(MATRIX + '/RNA.counts.matrix already exists, skip!')
    else:
        list_2_file(feature_counts_list, MATRIX +
                    '/feature_counts_list.txt')
        with ProcessPoolExecutor(max_workers=1) as executor:
            out_matrix = executor.submit(
                run_merge_featurecounts, 'RNA', MATRIX + '/feature_counts_list.txt')
            out_matrix = out_matrix.result()
            # 等待任务完成
            executor.shutdown(wait=True)
            # 移动前缀为'RNA'的文件到MATRIX_OUT_DIR
            for file in out_matrix:
                shutil.move(file, MATRIX)
    #### DESeq2 ####
    #### 之后再写 ####
    #### gatk remove duplicates ####
    rm_dup_threads = 8
    rm_dup_dict = {}
    mkdir(DUP_BAM)
    parallel = cacu_threads(threads, rm_dup_threads)
    with ProcessPoolExecutor(max_workers=parallel) as executor:
        for bam in samtools_sort_dict:
            rm_dup_out = bam + '_sorted_dup.bam'
            rm_dup_dict[bam] = rm_dup_out
            # 如果结果文件已经存在，跳过
            if os.path.exists(DUP_BAM + '/' + rm_dup_out + '.metrics.txt'):
                print(rm_dup_out + ' already exists, skip!')
                continue
            else:
                executor.submit(
                    gatk_remove_duplicates, samtools_sort_dict[bam], rm_dup_out)
        # 等待所有任务完成
        executor.shutdown(wait=True)
        for bam in rm_dup_dict:
            if os.path.exists(DUP_BAM +'/'+ rm_dup_dict[bam] + '.metrics.txt'):
                print(rm_dup_dict[bam] + ' already exists, skip!')
                continue
            else:
                shutil.move(rm_dup_dict[bam], DUP_BAM)
                shutil.move(rm_dup_dict[bam] + '.metrics.txt', DUP_BAM)
    index_threads = 1
    parallel = cacu_threads(threads, index_threads)
    with ProcessPoolExecutor(max_workers=parallel) as executor:
        for bam in rm_dup_dict:
            if os.path.exists(DUP_BAM +'/'+ bam + '_sorted_dup.bam.csi'):
                print(bam + '_sorted_dup.bam.csi' + ' already exists, skip!')
                continue
            else:
                executor.submit(samtools_index,DUP_BAM + '/' + rm_dup_dict[bam])
    #### gatk add read groups ####
    add_rg_threads = 8
    add_rg_dict = {}
    mkdir(ADD_HEADER)
    parallel = cacu_threads(threads, add_rg_threads)
    with ProcessPoolExecutor(max_workers=parallel) as executor:
        for bam in rm_dup_dict:
            add_rg_out = bam + '_sorted_dup_rg.bam'
            add_rg_dict[bam] = add_rg_out
            head_id = samples_dict[bam][2]
            # 如果结果文件已经存在，跳过
            if os.path.exists(VCF + '/all.vcf.gz'):
                print(add_rg_out + ' already exists, skip!')
                continue
            else:
                executor.submit(
                    gatk_add_read_groups, DUP_BAM + '/' + rm_dup_dict[bam], ADD_HEADER + '/' + add_rg_out, head_id)
    index_threads = 1
    parallel = cacu_threads(threads, index_threads)
    with ProcessPoolExecutor(max_workers=parallel) as executor:
        for bam in add_rg_dict:
            if os.path.exists(VCF + '/all.vcf.gz'):
                print(bam + '_sorted_dup_rg.bam.csi' + ' already exists, skip!')
                continue
            else:
                executor.submit(samtools_index,ADD_HEADER + '/' + add_rg_dict[bam])
    #### gatk split n trim ####
    split_n_trim_threads = 8
    gvcf_dict = {}
    mkdir(GVCF)
    parallel = cacu_threads(threads, split_n_trim_threads)
    with ProcessPoolExecutor(max_workers=parallel) as executor:
        for bam in add_rg_dict:
            bam_file = ADD_HEADER + '/' + add_rg_dict[bam]
            gvcf_dict[bam] = bam + '.g.vcf.gz'
            # 如果结果文件已经存在，跳过
            if os.path.exists(GVCF + '/' + bam + '.g.vcf.gz'):
                print(bam + '.g.vcf.gz' + ' already exists, skip!')
                continue
            else:
                executor.submit(
                    gatk_split, bam_file, genome, split_n_trim_threads,bam)
    gvcf_list = []
    for name in gvcf_dict:
        gvcf_list.append(GVCF + '/' + gvcf_dict[name])
        if not os.path.exists(GVCF + '/' + gvcf_dict[name]):
            shutil.move(gvcf_dict[name],GVCF)
            # 通配符匹配移动
            mv_file(gvcf_dict[name] + '.*' + 'g.vcf.gz', GVCF)
        else:
            print(GVCF + '/' + gvcf_dict[name] + ' already exists, skip!')

    index_threads = 1
    parallel = cacu_threads(threads, index_threads)
    with ProcessPoolExecutor(max_workers=parallel) as executor:
        for gvcf in gvcf_list:
            if os.path.exists(gvcf + '.csi') or os.path.exists(gvcf + '.tbi'):
                print(gvcf + '.csi' + ' already exists, skip!')
                continue
            else:
                executor.submit(bcfTools_index,gvcf)
    #### gatk merge gvcf ####
    merge_gvcf_threads = 8
    parallel = cacu_threads(threads, merge_gvcf_threads)
    with ProcessPoolExecutor(max_workers=parallel) as executor:
        if os.path.exists(GVCF + '/all.g.vcf.gz'):
            print(GVCF + '/all.g.vcf.gz' + ' already exists, skip!')
        else:
            executor.submit(gatk_merge_vcf,gvcf_list,genome,GVCF + '/all.g.vcf.gz')
    # bcfTools_index(GVCF + '/all.g.vcf.gz')
    if not os.path.exists(GVCF + '/all.g.vcf.gz' + '.csi'):
        if not os.path.exists(GVCF + '/all.g.vcf.gz' + '.tbi'):
            bcfTools_index(GVCF + '/all.g.vcf.gz')  
    #### gatk genotype gvcf to vcf ####
    mkdir(VCF)
    # genotype_gvcf_threads = 8
    # parallel = cacu_threads(threads, genotype_gvcf_threads)
    # with ProcessPoolExecutor(max_workers=parallel) as executor:
    #     if os.path.exists(VCF + '/all.vcf.gz'):
    #         print(VCF + '/all.vcf.gz' + ' already exists, skip!')
    #     else:
    #         executor.submit(gatk_genotype_vcf,GVCF + '/all.g.vcf.gz',genome,VCF + '/all.vcf.gz')
    #     # 等待任务完成
    #     executor.shutdown(wait=True)
    # 全部线程运行上面的任务
    if os.path.exists(VCF + '/all.vcf.gz'):
        print(VCF + '/all.vcf.gz' + ' already exists, skip!')
    else:
        gatk_genotype_vcf(GVCF + '/all.g.vcf.gz',VCF + '/all.vcf.gz',genome)
    
    if not os.path.exists(VCF + '/all.vcf.gz.csi'):
        if not os.path.exists(VCF + '/all.vcf.gz.tbi'):
            if os.path.exists(VCF + '/all.vcf.gz'):
                bcfTools_index(VCF + '/all.vcf.gz')
    bcf_calling_threads = 8
    bcf_list = []
    parallel = cacu_threads(threads, bcf_calling_threads)
    with ProcessPoolExecutor(max_workers=parallel) as executor:
        for bam in samtools_sort_dict:
            bam_file = samtools_sort_dict[bam]
            vcf_out = VCF + '/bcf_' + bam + '.vcf.gz'
            bcf_list.append(vcf_out)
            # 如果结果文件已经存在，跳过
            if os.path.exists('bcf_all.vcf.gz'):
                print(vcf_out + ' already exists, skip!')
                continue
            else:
                executor.submit(
                    bcfTools_call, bam_file, vcf_out,genome)
        # 等待任务完成
        executor.shutdown(wait=True)
    bcf_index_threads = 1
    parallel = cacu_threads(threads, bcf_index_threads)
    with ProcessPoolExecutor(max_workers=parallel) as executor:
        for vcf in bcf_list:
            if os.path.exists(vcf + '.csi') or os.path.exists(vcf + '.tbi'):
                print(vcf + '.csi' + ' already exists, skip!')
                continue
            else:
                executor.submit(bcfTools_index,vcf)
    if os.path.exists(VCF + '/bcf_all.vcf.gz'):
        print(VCF + '/bcf_all.vcf.gz' + ' already exists, skip!')
    else:
        bcfTools_merge(bcf_list, VCF + '/bcf_all.vcf.gz')
    if not os.path.exists(VCF + '/bcf_all.vcf.gz.csi'):
        if not os.path.exists(VCF + '/bcf_all.vcf.gz.tbi'):
            if os.path.exists(VCF + '/bcf_all.vcf.gz'):
                bcfTools_index(VCF + '/bcf_all.vcf.gz')
    # 删除dup bam fastp中间文件
    print('### remove dup bam fastp intermediate files ###')
    rm_files_list = []
    for sample in samtools_sort_dict:
        dup_out = DUP_BAM + '/' + sample + '_sorted_dup.bam'
        dupindex_out = DUP_BAM + '/' + sample + '_sorted_dup.bam.bai'
        sam_out = BAM_OUT_DIR + '/' + sample + '.sam'
        # bam_out = BAM_OUT_DIR + '/' + sample + '_sorted.bam'
        # bamindex_out = BAM_OUT_DIR + '/' + sample + '_sorted.bam.bai'
        bcf_vcf = VCF + '/bcf_' + sample + '.vcf.gz' 
        bcf_vcf_index = VCF + '/bcf_' + sample + '.vcf.gz.tbi'
        rm_files_list.append(bcf_vcf)
        rm_files_list.append(bcf_vcf_index)
        rm_files_list.append(dup_out)
        rm_files_list.append(dupindex_out)
        rm_files_list.append(sam_out)
        # rm_files_list.append(bam_out)
        # rm_files_list.append(bamindex_out)
    for fastq in fastp_dict:
        rm_files_list.append(fastp_dict[fastq][0])
        rm_files_list.append(fastp_dict[fastq][1])

    # 将bcf_all.vcf.gz 和 all.vcf.gz 转换为hapmap格式
    print('### convert bcf_all.vcf.gz and all.vcf.gz to hapmap format ###')
    if os.path.exists(VCF + '/bcf.hmp.txt'):
        print(VCF + '/bcf.hm.txt' + ' already exists, skip!')
    else:
        vcf_2_hapmap(VCF + '/bcf_all.vcf.gz',VCF + '/bcf')
    if os.path.exists(VCF + '/all.hmp.txt'):
        print(VCF + '/all.hmp.txt' + ' already exists, skip!')
    else:
        vcf_2_hapmap(VCF + '/all.vcf.gz',VCF + '/all')
    rm_files(rm_files_list)
    print('all done!')
    #### gatk variant filter ####
    # print("### gatk remove duplicates ###")
    # mkdir(BAM_OUT_DIR)
    # gatk_remove_duplicates_threads = 8
    # dup_dict = {}
    # with ProcessPoolExecutor(max_workers=gatk_remove_duplicates_threads) as executor:
    #     for sample in samtools_sort_dict:
    #         dup_out = BAM_OUT_DIR + '/' + sample + '_sorted_dup.bam'
    #         dup_dict[sample] = dup_out
    #         # 如果结果文件已经存在，跳过
    #         if os.path.exists(dup_out):
    #             print(dup_out + ' already exists, skip!')
    #             continue
    #         else:
    #             executor.submit(gatk_remove_duplicates,
    #                             samtools_sort_dict[sample], dup_out)



if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me! ;-) See you!\n")
        sys.exit(0)
