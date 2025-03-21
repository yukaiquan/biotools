#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse


def samtools_view(chr, reference, start, end, output):
    cmd = f"samtools faidx {reference} {chr}:{start}-{end} > {output}"
    return cmd

# 执行cmd命令
def execute_cmd(cmd):
    import subprocess
    subprocess.run(cmd, shell=True, check=True)
    return True

def main():
    # 输入参数
    parser = argparse.ArgumentParser(description='samtools faidx please input reference and samples file, the format of samples file is chr\tstart\tend')
    parser.add_argument('-r', '--reference', type=str, required=True, help='reference file')
    parser.add_argument('-s','--samples', type=str, required=True, help='samples file')
    args = parser.parse_args()
    reference = args.reference
    samples = args.samples
    print("yukaiquan 1962568272@qq.com")
    print(f"reference: {reference}")
    print(f"samples: {samples}")

    # 读取样本文件
    cmd_list = []
    with open(samples, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            # 删除空行
            if line.strip() == '':
                continue
            lines = line.strip().split()
            if len(lines) != 3:
                print(f"Error: invalid line {line}")
                exit(1)
            # print(lines)
            chr_name = lines[0]
            start = lines[1]
            end = lines[2]
            if int(start) > int(end):
                print(f"Error: start {start} is greater than end {end}")
                exit(1)
            output = chr_name + '_' + start + '_' + end + '.fa'
            cmd = samtools_view(chr_name, reference, start, end, output)
            print(cmd)
            cmd_list.append(cmd)

    # 执行命令
    for cmd in cmd_list:
        execute_cmd(cmd)
    
    print("Done!")
if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        print(f"Error: {e}")
        exit(1)

