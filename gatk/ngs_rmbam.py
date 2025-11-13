#!/usr/bin/env python3
import sys
import os
from subprocess import Popen, PIPE
from concurrent.futures import ProcessPoolExecutor
import argparse
# samtools 版本: 1.12 太低了不支持bwa mem的输出顺便排序


# GATK_SIF = 'singularity exec ~/soft/gatk_461.sif gatk'
GATK_SIF = 'gatk'
SAMTOOLS = '/work/home/acfaa2ssz7/soft/samtools-1.10/samtools'
DEPTH = '/work/home/acfaa2ssz7/soft/pandepth'





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



class BAM:
    def __init__(self, input_bam: str, output_bam: str, thread: int):
        self.input_bam = input_bam
        self.output_bam = output_bam
        self.thread = thread

    def rmdup(self) -> str:
        # gatk MarkDuplicates -I 02wumang.bam -O 02wumang.rmdup.bam -M 02wumang.rmdup.metrics.txt
        gatk_markdup = GATK_SIF + ' MarkDuplicates ' + ' -I ' + self.input_bam + ' -O ' + self.output_bam + ' -M ' + self.output_bam + '.metrics.txt'
        return gatk_markdup
    def samtools_index(self) -> str:
        # 判断
        if not os.path.exists(self.output_bam):
            return 'NA'
        # samtools index 单线程
        samtools_index = SAMTOOLS + ' index -@ ' + str(self.thread) + ' -c ' + self.output_bam
        return samtools_index
    def depth_of_coverage(self) -> str:
        # if not os.path.exists(self.output_bam + '.chr.stat.gz'):
        #     return 'NA'
        # depth_of_coverage
        depth = DEPTH  + ' -i ' + self.output_bam + ' -t ' + str(self.thread) + ' -o ' + self.output_bam
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
    parser.add_argument('-i', '--input', type=str, required=True, help='input bam file')
    parser.add_argument('-o', '--output', type=str, required=True, help='output rmdup bam file')
    parser.add_argument('-t', '--thread', type=int, default=4, help='thread number(default: 4)')
    parser.add_argument('-d', '--depth', type=bool, default=False, help='calculate depth of coverage(default: False)')
    args = parser.parse_args()

    print('#' * 60 + '\n' + ' ' * 20 + 'bam rmdup' + '\n' + '#' * 60 + '\n')
    # check file

    output_bam = args.output
    input_bam = args.input
    thread = args.thread
    depth = args.depth
    check_file(input_bam)
    bam = BAM(input_bam, output_bam, thread)

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
    # 检查输出文件是否存在，存在则跳过

    if not os.path.exists(output_bam):
        rm_bam_dup = bam.rmdup()
        print(rm_bam_dup)
        command = Command(rm_bam_dup, args.thread)
        out_log = command.run_with_return()
        print(out_log)
    # 检查输出文件是否大小大于0
    if os.path.exists(output_bam) and os.path.getsize(output_bam) > 0:
        print('Output file ' + output_bam + ' is ok!')
    else:
        print('Output file ' + output_bam + ' is empty or not exists!')
        sys.exit(1)
    # 
    ## 单线程
    ## # samtools index
    ## samtools_index = 'samtools index ' + args.output + '/' + args.prefix + '.rmdup.bam'
    ## print(samtools_index)
    ## os.system(samtools_index)
    # 检查*.rmdup.bam.csi文件是否存在
    if not os.path.exists(output_bam + '.csi'):
        samtools_rm_dup_index = bam.samtools_index()
        print(samtools_rm_dup_index)
        command = Command(samtools_rm_dup_index, args.thread)
        out_log = command.run_with_return()
        print(out_log)
    else:
        print('Output file ' + output_bam + '.csi is ok!')

    # 检查*.rmdup.bam.csi文件是否大小大于0
    if os.path.exists(output_bam + '.csi') and os.path.getsize(output_bam + '.csi') > 0:
        print('Output file ' + output_bam + '.csi is ok!')
        rm_file(input_bam)
        print('Input file ' + input_bam + ' is removed!')
    else:
        print('Output file ' + output_bam + '.csi is empty or not exists!')
        sys.exit(1)


    if not depth:
        print('No need to calculate depth of coverage!')
        print('#' * 60 + '\n' + ' ' * 20 + 'bam rmdup done!' + '\n' + '#' * 60 + '\n')
        return None
    else:
        depth_shell = bam.depth_of_coverage()
        print(depth_shell)
        command = Command(depth_shell, args.thread)
        out_log = command.run_with_return()
        print(out_log)



        # rm_file(args.output + '/' + args.prefix + '.bam')
        print('#' * 60 + '\n' + ' ' * 20 + 'bam rmdup done!' + '\n' + '#' * 60 + '\n')
        return None

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write('User interrupt me! ;-) See you!\n')
        sys.exit(0)
