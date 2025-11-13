#!/usr/bin/env python3
import sys
import os
from subprocess import Popen, PIPE
from concurrent.futures import ProcessPoolExecutor
import argparse
# samtools 版本: 1.12 太低了不支持bwa mem的输出顺便排序

# GATK_SIF = 'singularity exec /public/share/acfaa2ssz7/soft/singularity/gatk45.sif gatk'
# GATK_SIF = 'singularity exec /public/home/pengyuanying/soft/gatk_461.sif gatk'

# GATK_SIF = '/public/home/pengyuanying/.conda/envs/bwa/bin/gatk'
GATK_SIF = '/work/share/acfaa2ssz7/miniconda3/envs/gatk4/bin/gatk'
BWATOOLS = '/work/share/acfaa2ssz7/miniconda3/envs/bwa/bin/bwa'
# BWATOOLS = '/public/home/pengyuanying/.conda/envs/bwa/bin/bwa'
SAMTOOLS = '/work/home/acfaa2ssz7/soft/samtools-1.10/samtools'
# SAMTOOLS = '/public/home/pengyuanying/.conda/envs/bwa/bin/samtools'
FASTP = '/work/share/acfaa2ssz7/miniconda3/envs/bwa/bin/fastp'
# FASTP = '/public/home/pengyuanying/.conda/envs/bwa/bin/fastp'
DEPTH = '/work/home/acfaa2ssz7/soft/pandepth'
# DEPTH = '/public/home/pengyuanying/soft/pandepth'




def check_file(file_path: str) -> None:
    if not os.path.exists(file_path):
        print('Error: file or directory not exists ' + file_path)
        sys.exit(1)
    return None

def check_dir(dir_path: str) -> None:
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
    return None

def rm_file(file_path: str) -> None:
    if os.path.exists(file_path):
        os.remove(file_path)
    return None

class Fastp:
    def __init__(self, input1: str, input2: str, output_dir: str, prefix: str, thread: int):
        self.input1 = input1
        self.input2 = input2
        self.output_dir = output_dir
        self.prefix = prefix
        self.thread = thread if thread <= 16 else 16
        self.clean_1 = f"{output_dir}/{prefix}_clean_1.fq.gz"
        self.clean_2 = f"{output_dir}/{prefix}_clean_2.fq.gz"
        self.html_report = f"{output_dir}/{prefix}_fastp.html"
        self.json_report = f"{output_dir}/{prefix}_fastp.json"

    def build_command(self) -> str:
        # 如果存在self.html_report就返还'NA'
        if os.path.exists(self.html_report):
            return 'NA'
        # 构建fastp命令
        cmd = (
            f"{FASTP} -i {self.input1} -I {self.input2} "
            f"-o {self.clean_1} -O {self.clean_2} "
            f"--html {self.html_report} --json {self.json_report} "
            f"--thread {self.thread}"
        )
        return cmd

class BWA:
    def __init__(self, reference: str, thread: int, flag: str, output: str, prefix: str):
        self.reference = reference
        self.thread = thread
        self.flag = flag
        self.output = output
        self.prefix = prefix

    def bwa_mem(self, input1: str, input2: str) -> str:
        # bwa mem
        # bwa mem /mnt/disk4/sicau/BN220729WH01S05N5_3Samples/new_genome/ACD/bwaindex/ACD.sativa.fasta 0002.out.wumang_1.fq.gz 0002.out.wumang_2.fq.gz -t 16 -R "@RG\tID:wumang\tLB:wumang\tPL:ILLUMINA\tSM:wumang" -o 02wumang.sam 
        bwa_flag = '@RG\\tID:' + self.flag + '\\tLB:' + self.flag + '\\tPL:ILLUMINA\\tSM:' + self.flag
        bwa_mem = BWATOOLS + ' mem ' + self.reference + ' ' + input1 + ' ' + input2 + ' -K 1000000000 -t ' + str(self.thread) + ' -R "' + bwa_flag + '" -o ' + self.output + '/' + self.prefix + '.sam'
        return bwa_mem

    def samtools_sort(self) -> str:
        # samtools sort
        samtools_sort = SAMTOOLS + ' sort -@ ' + str(self.thread)  + ' -o ' + self.output + '/' + self.prefix + '.bam ' + self.output + '/' + self.prefix + '.sam'
        return samtools_sort
    
    def bwa_mem_pipeline(self, input1: str, input2: str) -> str:
        # 整合bwa比对+samtools转换+samtools排序的完整管道命令
        # 使用双引号包裹RG参数避免转义问题
        bwa_flag = f"@RG\\tID:{self.flag}\\tLB:{self.flag}\\tPL:ILLUMINA\\tSM:{self.flag}"
        
        command = (
            f"{BWATOOLS} mem -K 1000000000 -t {self.thread} "  # BWA参数
            f"-R \"{bwa_flag}\" "  # Read Group信息
            f"{self.reference} {input1} {input2} "  # 输入文件
            "| "  # 管道传递
            f"{SAMTOOLS} view -@ 4 -Sb "  # 转换BAM
            "| "  # 继续管道传递
            f"{SAMTOOLS} sort -@ 4 -m 4G "  # 排序
            f"-o {self.output}/{self.prefix}.bam "  # 输出路径
            # "-"  # 从标准输入读取
        )
        return command

    def samtools_index(self) -> str:
        # samtools index 单线程
        samtools_index = SAMTOOLS + ' index -@ ' + str(self.thread) + ' -c ' + self.output + '/' + self.prefix + '.bam'
        return samtools_index

    def samtools_flagstat(self) -> str:
        # samtools flagstat 单线程
        samtools_flagstat = SAMTOOLS + ' flagstat -@ ' + str(self.thread) + ' ' + self.output + '/' + self.prefix + '.bam > ' + self.output + '/' + self.prefix + '.bam.flagstat'
        return samtools_flagstat

    def rm_dup(self) -> str:
        # 低线程占内存
        # rm_dup
        # gatk MarkDuplicates -I 02wumang.bam -O 02wumang.rmdup.bam -M 02wumang.rmdup.metrics.txt
        gatk_markdup = GATK_SIF + ' MarkDuplicates ' + ' -I ' + self.output + '/' + self.prefix + '.bam -O ' + self.output + '/' + self.prefix + '.rmdup.bam -M ' + self.output + '/' + self.prefix + '.rmdup.metrics.txt'
        return gatk_markdup
    def samtools_index_rmdup(self) -> str:
        # 判断
        if not os.path.exists(self.output + '/' + self.prefix + '.rmdup.bam'):
            return 'NA'
        # samtools index 单线程
        samtools_index = SAMTOOLS + ' index -@ ' + str(self.thread) + ' -c '  + self.output + '/' + self.prefix + '.rmdup.bam'
        return samtools_index
    
    def depth_of_coverage(self) -> str:
        # if not os.path.exists(self.output + '/' + self.prefix + '.chr.stat.gz'):
        #     return 'NA'
        # depth_of_coverage
        depth = DEPTH  + ' -i ' +self.output + '/' + self.prefix + '.rmdup.bam' +  ' -t ' + str(self.thread) + ' -o ' + self.output + '/' + self.prefix
        return depth

class Command:
    def __init__(self, command: str, all_thread: int):
        self.command = command
        self.all_thread = all_thread

    # 私有方法
    def __caculate_parallel_thread(self,thread: int) -> int: 
        if self.all_thread <= thread:
            parallel = 1
            thread = self.all_thread
        elif self.all_thread > thread and self.all_thread % thread == 0:
            parallel = self.all_thread // thread
        else:
            parallel = self.all_thread // thread + 1
        return parallel

    # 普通命令
    def run(self) -> None:
        # print(self.command)
        os.system(self.command)
        return None
    # 带返回值的命令
    def run_with_return(self) -> str:
        # print(self.command)
        # return os.popen(self.command).read()
        # 等待命令执行完成
        if self.command == 'NA':
            return 'Fastp pass !'
        process = Popen(self.command, shell=True, stdout=PIPE, stderr=PIPE)
        stdout, stderr = process.communicate()
        print(stdout.decode('utf-8'))
        print(stderr.decode('utf-8'))
        process.wait()
        return stdout.decode('utf-8')
    # 并行命令
    def run_parallel(self,thread: int) -> str:
        print(self.command)
        parallel = self.__caculate_parallel_thread(thread)
        with ProcessPoolExecutor(max_workers=parallel) as executor:
            process = executor.submit(os.popen, self.command)
            return process.result().read()

def main():
    parser = argparse.ArgumentParser(description='bwa mem')
    parser.add_argument('-i', '--input', type=str, required=True, help='input fastq first file')
    parser.add_argument('-j', '--input2', type=str, required=True, help='input fastq second file')
    parser.add_argument('-o', '--output', type=str, required=True, help='output sam file')
    parser.add_argument('-r', '--reference', type=str, required=True, help='reference genome')
    parser.add_argument('-t', '--thread', type=int, default=8, help='thread number')
    parser.add_argument('-f', '--flag', type=str, help='bwa flag')
    parser.add_argument('-p', '--prefix', type=str, default='bwa', help='output prefix')
    args = parser.parse_args()

    print('#' * 60 + '\n' + ' ' * 20 + 'bwa mem' + '\n' + '#' * 60 + '\n')
    # check file
    check_file(args.input)
    check_file(args.input2)
    check_file(args.reference)
    check_dir(args.output)
        # 创建Fastp实例并运行fastp
    fastp_instance = Fastp(
        input1=args.input,
        input2=args.input2,
        output_dir=args.output,
        prefix=args.prefix,
        thread=args.thread
    )
    print(fastp_instance.build_command())
    command = Command(fastp_instance.build_command(), args.thread)
    out_log = command.run_with_return()
    print(out_log)
    clean_1_fastq = args.output + '/' + args.prefix + "_clean_1.fq.gz"
    clean_2_fastq = args.output + '/' + args.prefix + "_clean_2.fq.gz"
    # bwa mem
    # bwa mem /mnt/disk4/sicau/BN220729WH01S05N5_3Samples/new_genome/ACD/bwaindex/ACD.sativa.fasta 0002.out.wumang_1.fq.gz 0002.out.wumang_2.fq.gz -t 16 -R "@RG\tID:wumang\tLB:wumang\tPL:ILLUMINA\tSM:wumang" -o 02wumang.sam 
    # bwa_flag = '@RG\\tID:' + args.flag + '\\tLB:' + args.flag + '\\tPL:ILLUMINA\\tSM:' + args.flag
    # bwa_mem = 'bwa mem ' + args.reference + ' ' + args.input + ' ' + args.input2 + ' -t ' + str(args.thread) + ' -R "' + bwa_flag + '" -o ' + args.output + '/' + args.prefix + '.sam'
    # print(bwa_mem)
    # os.system(bwa_mem)
    # bwa_thread = 32
    bwa_mem = BWA(args.reference, args.thread, args.flag,args.output, args.prefix)
    # # 生成命令
    # cmd = bwa_mem.bwa_mem(clean_1_fastq, clean_2_fastq)
    # print(cmd)
    # # 执行命令
    # command = Command(cmd, args.thread)
    # out_log = command.run_with_return()
    # print(out_log)

    # # samtools sort
    # # samtools_sort = 'samtools sort -@ ' + str(args.thread)  + ' -o ' + args.output + '/' + args.prefix + '.bam ' + args.output + '/' + args.prefix + '.sam'
    # # print(samtools_sort)
    # # os.system(samtools_sort)
    # samtools_sort = bwa_mem.samtools_sort()
    # print(samtools_sort)
    # command = Command(samtools_sort, args.thread)
    # out_log = command.run_with_return()
    
    cmd = bwa_mem.bwa_mem_pipeline(clean_1_fastq, clean_2_fastq)
    print(cmd)
    command = Command(cmd, args.thread)
    out_log = command.run_with_return()
    print(out_log)

    # 检查输出的bam文件是否大于1KB
    output_bam = args.output
    if not os.path.exists(output_bam) or os.path.getsize(output_bam) < 1024:
        print('Error: Output BAM file is empty or not exists: ' + output_bam)
        sys.exit(1)
    else:
        print('Output BAM file is ok: ' + output_bam)
        # 删除fastq文件
        # rm_file(clean_1_fastq)
        # rm_file(clean_2_fastq)
        rm_file(args.input)
        rm_file(args.input2)
        print('Cleaned up input fastq files.')

    # 单线程
    # samtools index
    # samtools_index = 'samtools index ' + args.output + '/' + args.prefix + '.bam'
    # print(samtools_index)
    # os.system(samtools_index)
    #samtools_index = bwa_mem.samtools_index()
    #print(samtools_index)
    #command = Command(samtools_index, args.thread)
    #out_log = command.run_with_return()
    #print(out_log)

    ## 单线程
    ## # samtools flagstat
    ## samtools_flagstat = 'samtools flagstat ' + args.output + '/' + args.prefix + '.bam > ' + args.output + '/' + args.prefix + '.bam.flagstat'
    ## print(samtools_flagstat)
    ## os.system(samtools_flagstat)
    #samtools_flagstat = bwa_mem.samtools_flagstat()
    #print(samtools_flagstat)
    #command = Command(samtools_flagstat, args.thread)
    #out_log = command.run_with_return()
    #print(out_log)
    ## rm sam
    ## rm_file(args.output + '/' + args.prefix + '.sam')
    ## 低线程占内存
    ## rm_dup
    ## rm_dup = rm_bam_dup(args.output + '/' + args.prefix + '.bam', args.thread, args.output + '/' + args.prefix + '.rmdup.bam')
    ## print(rm_dup)
    ## os.system(rm_dup)
    #rm_bam_dup = bwa_mem.rm_dup()
    #print(rm_bam_dup)
    #command = Command(rm_bam_dup, args.thread)
    #out_log = command.run_with_return()
    #print(out_log)

    ## 单线程
    ## # samtools index
    ## samtools_index = 'samtools index ' + args.output + '/' + args.prefix + '.rmdup.bam'
    ## print(samtools_index)
    ## os.system(samtools_index)
    #samtools_rm_dup_index = bwa_mem.samtools_index_rmdup()
    #print(samtools_rm_dup_index)
    #command = Command(samtools_rm_dup_index, args.thread)
    #out_log = command.run_with_return()
    #print(out_log)

    ## rm_file(args.output + '/' + args.prefix + '.bam')
    #depth_shell = bwa_mem.depth_of_coverage()
    #print(depth_shell)
    #command = Command(depth_shell, args.thread)
    #out_log = command.run_with_return()
    #print(out_log)




    print('#' * 60 + '\n')
    print('BWA MEM DONE!\n')
    print('#' * 60 + '\n')
    return None

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write('User interrupt me! ;-) See you!\n')
        sys.exit(0)
