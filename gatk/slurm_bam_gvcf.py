#!/usr/bin/env python3
#-- coding:UTF-8 --
import sys
import os
import time
import glob
from subprocess import Popen, PIPE
from concurrent.futures import ProcessPoolExecutor
import argparse
import concurrent
import subprocess

# GATK: str = '/public/home/acfaa2ssz7/soft/gatk4.py'
# GATK_SIF = 'singularity exec /public/share/acfaa2ssz7/soft/singularity/gatk45.sif gatk'
# GATK_SIF = 'singularity exec /work/home/acfaa2ssz7/soft/gatk_461.sif gatk'
GATK_SIF = 'gatk'

CHR_LIST = ['chr1A', 'chr2A', 'chr3A', 'chr4A', 'chr5A', 'chr6A', 'chr7A',
                  'chr1C', 'chr2C', 'chr3C', 'chr4C', 'chr5C', 'chr6C', 'chr7C',
                  'chr1D', 'chr2D', 'chr3D', 'chr4D', 'chr5D', 'chr6D', 'chr7D', 'chrUn']
# SAMTOOLS = '/public/home/pengyuanying/.conda/envs/bwa/bin/samtools'
SAMTOOLS = '/work/home/acfaa2ssz7/soft/samtools-1.10/samtools'
BCFTOOLS = '/work/share/acfaa2ssz7/miniconda3/envs/bwa/bin/bcftools'
# BCFTOOLS = '/public/home/pengyuanying/.conda/envs/bwa/bin/bcftools'

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



# def gatk4(ref: str, input_bam: str, output: str, threads: int) -> int:
#     cmd = f'python {GATK} -P {threads} -I {input_bam} -R {ref} -O {output}'
#     p = Popen(cmd, shell=True)
#     return_code = p.wait()  # 等待子进程结束，并返回状态码；0表示正常结束，非0表示异常结束
#     return return_code


def gatk4_merge(chrom_gvcf, output):
    cmd = f'{GATK_SIF} GatherVcfs -RI TRUE -I {" -I ".join(chrom_gvcf)} -O {output}'
    p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    return_code = p.wait()
    return return_code


def samtools_index(bam: str) -> True:
    '''
    samtools index -c bam.file
    '''
    if os.path.exists(f'{bam}.csi'):
        return True
    cmd = f'{SAMTOOLS} index -@ 4 -c {bam}'
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
    # cmd = f'bcftools index -t {gvcf}'
    cmd = f'{BCFTOOLS} index -c {gvcf}'
    print(cmd)
    # if os.path.exists(f'{gvcf}.tbi'):
    #     return True
    if os.path.exists(f'{gvcf}.csi'):
        return True
    p = Popen(cmd, shell=True)
    return_code = p.wait()
    return return_code


def geno_type_GVCFs(ref: str, input_gvcf: str, output_vcf: str) -> bool:
    cmd = f'{GATK_SIF} GenotypeGVCFs -OVI False -R {ref} -V {input_gvcf} -O {output_vcf}'
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

def generate_sbatch_file(patient_id, job_name, job_out, job_err, job_mem, job_N, job_cpu, job_script,folder)->str:
    """
    Generate a sbatch file for submitting a job to the cluster.
    """
    
    sbatch_file_mem = f"""#!/bin/bash
#SBATCH -J {job_name}
#SBATCH -o {job_out}
#SBATCH -e {job_err}
#SBATCH -p {patient_id}
#SBATCH -N {job_N}
#SBATCH -n {job_cpu}
#SBATCH --mem={job_mem}

module purge

source /public/home/acfaa2ssz7/.bashrc
conda activate gatk4

{job_script}
wait
"""
    sbatch_file = f"""#!/bin/bash
#SBATCH -J {job_name}
#SBATCH -o {job_out}
#SBATCH -e {job_err}
#SBATCH -p {patient_id}
#SBATCH -N {job_N}
#SBATCH -n {job_cpu}

module purge

source /public/home/acfaa2ssz7/.bashrc
conda activate gatk4

{job_script}
wait
"""
    if job_mem == "0":
        with open(f"{folder}/{job_name}.sbatch", "w") as f:
            f.write(sbatch_file)
    else:
        with open(f"{folder}/{job_name}.sbatch", "w") as f:
            f.write(sbatch_file_mem)
    return f"{folder}/{job_name}.sbatch"

def genetion_gvcf(ref: str, input_bam: str, output_gvcf: str, chrom: str, process:str,job:str,node:str,threads:int) -> str:
    out_put_name = output_gvcf + '.' + chrom + '.g.vcf.gz'
    # 判断文件是否已经存在
    if os.path.exists(out_put_name):
        return "NA"
    else:
        # if node == 0:
        #     cmd = f'srun -J {job}_{chrom} -p {process} -n {str(threads)} --mem=2000MB -o std.out.%j -e std.err.%j {GATK_SIF} --java-options -Xmx6G HaplotypeCaller -I {input_bam} -O {out_put_name} -R {ref} --emit-ref-confidence GVCF -OVI False -L {chrom} &'
        # else:
        #     cmd = f'srun -J {job}_{chrom} -p {process} -N {str(node)} -n {str(threads)} --mem=2000MB -o std.out.%j -e std.err.%j {GATK_SIF} --java-options -Xmx6G HaplotypeCaller -I {input_bam} -O {out_put_name} -R {ref} --emit-ref-confidence GVCF -OVI False -L {chrom} &'
        cmd = f'{GATK_SIF} --java-options -Xmx6G HaplotypeCaller -I {input_bam} -O {out_put_name} -R {ref} --emit-ref-confidence GVCF -OVI False -L {chrom}'
        print(cmd)
        return cmd


def gvcf_cmd_list(ref, input_bam,folder, output_gvcf, process:str,job:str,node:str,threads:int,mem:str) -> dict:
    # sbatch_file_list = []
    sbatch_file_dict = {}
    for chrom in CHR_LIST:
        cmd = genetion_gvcf(ref, input_bam, output_gvcf, chrom, process,job,node,threads)
        gvcf_file = output_gvcf + '.' + chrom + '.g.vcf.gz'
        if cmd == "NA":
            continue
        else:
            job_name = job + '_' + chrom
            job_out = f'{output_gvcf}_{job_name}.out.%j'
            job_err = f'{output_gvcf}_{job_name}.err.%j'
            job_mem = mem
            job_N = node
            job_cpu = threads
            job_script = cmd
            # 生成sbatch文件
            sbatch_file = generate_sbatch_file(process, job_name, job_out, job_err, job_mem, job_N, job_cpu, job_script,folder)
            # sbatch_file_list.append(sbatch_file)
            if gvcf_file not in sbatch_file_dict:
                sbatch_file_dict[gvcf_file] = sbatch_file
            else:
                print(f'Error: {gvcf_file} already exists in {sbatch_file_dict[gvcf_file]}')
                sys.exit(1)
    return sbatch_file_dict


def pysh(cmd) -> bool:
    '''
    执行输入的shell命令,并返回执行状态bool
    '''
    # p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    # stdout, stderr = p.communicate()
    # return_code = p.wait()
    if "HaplotypeCaller" in cmd:
        output_gvcf = cmd.split('-O ')[1].split('-R')[0].strip()
        if os.path.exists(output_gvcf):
            return True
    subprocess.run(cmd, shell=True)
    # 直接运行不等待任务返回
    # process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # # 获取子进程的输出
    # output, error = process.communicate()

    # # 打印输出
    # print("输出：", output.decode())
    # print("错误：", error.decode())
    # if return_code != 0:
    #     print(f'Error: {cmd} {stderr}')
    #     sys.exit(1)
    return True


def gatk4_merge(chrom_gvcf: list, output: str) -> int:
    cmd = f'{GATK_SIF} GatherVcfs -RI TRUE -I {" -I ".join(chrom_gvcf)} -O {output}'
    print(cmd)
    p = Popen(cmd, shell=True)
    return_code = p.wait()
    return return_code


def sample_fasq_list(fastq1_list: list, fastq2_list: list, sample_input_fq: str) -> list:
    with open(sample_input_fq, 'r') as f:
        for line in f:
            fastq1_list.append(line.strip().split()[0])
            fastq2_list.append(line.strip().split()[1])
    return fastq1_list, fastq2_list


def file_exists(file: str) -> bool:
    if os.path.exists(file):
        return True
    else:
        return False

def read_bam_samples(bam_samples_list: str) -> dict:
    samples_dict = {}
    with open(bam_samples_list, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            line_list = line.strip().split(',')
            sample_name = line_list[0]
            # 判断文件是否存在
            if file_exists(line_list[1]):
                samples_dict[sample_name] = line_list[1]
            else:
                print(f'Error: {line_list[1]} not exists')
                sys.exit(1)
    return samples_dict

def read_fai(fai_file: str) -> list:
    chr_length_dict = {}
    if file_exists(fai_file):
        with open(fai_file, 'r') as f:
            for line in f:
                line_list = line.strip().split()
                chr_name = line_list[0]
                chr_length = line_list[1]
                if chr_name in chr_length_dict:
                    print(f'Error: {chr_name} exists')
                    sys.exit(1)
                else:
                    chr_length_dict[chr_name] = chr_length
    else:
        print(f'Error: {fai_file} not exists')
        sys.exit(1)
    # 根据染色体长度排序
    chr_tup = sorted(chr_length_dict.items(), key=lambda x: int(x[1]), reverse=True)
    chr_list = [i[0] for i in chr_tup]
    return chr_list

def sorted_chr(chr_len_sorted:list,cmd_list:list)->list:
    # 为了充分运用计算资源，我们倾向于先计算较长染色体
    # print(chr_len_sorted)
    new_cmd_list = []
    for i in range(len(chr_len_sorted)):
        chr_name = chr_len_sorted[i]
        for j in range(len(cmd_list)):
            if chr_name in cmd_list[j]:
                new_cmd_list.append(cmd_list[j])
    return new_cmd_list

def main():
    parser = argparse.ArgumentParser(description='bam 2 vcf')
    parser.add_argument('-i', '--input', type=str, required=True, help='input rm dup bam')
    parser.add_argument('-r', '--reference', type=str, required=True, help='reference genome')
    parser.add_argument('-t', '--thread', type=int, default=2, help='thread number')
    parser.add_argument('-p','--procession', type=str,required=True, help='procession name')
    parser.add_argument('-m', '--mem', type=str, default='2G', help='memory')
    parser.add_argument('-j', '--job', type=str,required=True, help='job name')
    parser.add_argument('-n','--node', type=int,required=True, help='node numbers')
    parser.add_argument('-d', '--debug', type=bool, default=False, help='debug')
    parser.add_argument('-o', '--output', type=str,required=True, help='output vcf name')
    args = parser.parse_args()
    input_bam = args.input
    reference = args.reference
    thread = args.thread
    output_name = args.output
    process = args.procession
    job = args.job
    node = args.node
    mem = args.mem
    debug = args.debug
    # 判断文件是否存在
    if file_exists(input_bam):
        pass
    else:
        print(f'Error: {input_bam} not exists')
        sys.exit(1)
    # 判断文件是否存在
    if file_exists(reference):
        pass
    else:
        print(f'Error: {reference} not exists')
        sys.exit(1)
    folder = os.path.dirname(input_bam)
    output_path_name = os.path.join(folder, output_name)
    gvcf_chr = {}
    output_name_list = [output_path_name + '.' + chr +
                            '.g.vcf.gz' for chr in CHR_LIST]
    gvcf_chr[output_path_name] = output_name_list
    # 总的gvcf文件
    gvcf_file = output_path_name + '.g.vcf.gz'
    # 判断是否存在如果存在则退出
    if file_exists(gvcf_file):
        print(f'Error: {gvcf_file} exists')
        sys.exit(1)
    # 判断文件是否存在
    sbath_file_dict = gvcf_cmd_list(reference, input_bam,folder, output_path_name,process,job,node,thread,mem)
    # 直接运行
    for gvcf_name in sbath_file_dict.keys():
        sbatch_file = sbath_file_dict[gvcf_name]
        # 判断结果文件是否存在，如果存在则不
        if file_exists(gvcf_name):
            print(f'Error: {gvcf_name} exists')
            continue
        else:
            if debug == True:
                continue
            time.sleep(5)
            cmd = f'sbatch {sbatch_file}'
            res = pysh(cmd)
            print(f'Submitted job for tag: {gvcf_name}')
    if debug == True:
        pass
    else:
        print("All jobs submitted successfully.")
    
        print(f'{output_name} is generating; next you will merge gvcf files!')
    

if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        print(e)
        sys.exit(1)
