#!/usr/bin/env python
import sys
import os
import time
import glob
from subprocess import Popen, PIPE
from concurrent.futures import ProcessPoolExecutor
import argparse


def parse_args():
    parser = argparse.ArgumentParser(description='diamond blast.')
    parser.add_argument('-i', '--input', type=str,
                        required=True, help='Input genome pep 1')
    parser.add_argument('-I', '--input2', type=str,
                        required=True, help='Input genome pep 2')
    parser.add_argument('-o', '--output', type=str,
                        required=True, help='Output folder')
    parser.add_argument('-t', '--threads', type=int, default=16,
                        help='Number of threads')
    parser.add_argument('-s', '--self1', type=bool, default=False,
                        help='self blast 1')
    parser.add_argument('-S', '--self2', type=bool, default=False,
                        help='self blast 2')
    return parser.parse_args()


def mkdir(path):
    folder = os.path.exists(path)
    if not folder:  # 判断是否存在文件夹如果不存在则创建为文件夹
        os.makedirs(path)  # makedirs 创建文件时如果路径不存在会创建这个路径


def makeblastdb(input_file, output_file) -> int:
    """
    makeblastdb
    """
    cmd = 'diamond makedb --in {} -d {}'.format(input_file, output_file)
    # 检查input_file.dmnd是否存在
    if os.path.exists(output_file + '.dmnd'):
        return 0
    p = Popen(cmd, shell=True)
    return_code = p.wait()
    print(cmd)
    if return_code != 0:
        print(f'Error: {cmd}')
        sys.exit(1)
    return return_code


def blastp_run(input_file, output_file, db_file, threads) -> int:
    """
    blastp
    """
    cmd = f'diamond blastp -q {input_file} -d {db_file} -o {output_file} --more-sensitive -p {threads} --quiet -e 0.001'
    p = Popen(cmd, shell=True)
    return_code = p.wait()
    print(cmd)
    if return_code != 0:
        print(f'Error: {cmd}')
        sys.exit(1)
    return return_code


def main():
    """
    main
    """
    start_time = time.time()
    args = parse_args()
    output_dir = args.output
    input = args.input
    input2 = args.input2
    threads = args.threads
    self1 = args.self1
    self2 = args.self2
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    if makeblastdb(input, input) != 0:
        print('Error: makeblastdb genome 1')
        sys.exit(1)
    if makeblastdb(input2, input2) != 0:
        print('Error: makeblastdb genome 2')
        sys.exit(1)
    blastp_list = []
    blastp_list.append([input, output_dir + '/' + input +
                       '_' + input2 + '.blast', input2, threads])
    if self1:
        blastp_list.append([input, output_dir + '/' + input +
                            '_' + input + '.blast', input, threads])
    blastp_list.append([input2, output_dir + '/' + input2 +
                       '_' + input + '.blast', input, threads])
    if self2:
        blastp_list.append([input2, output_dir + '/' + input2 +
                            '_' + input2 + '.blast', input2, threads])
    for blastp in blastp_list:
        if blastp_run(blastp[0], blastp[1], blastp[2], blastp[3]) != 0:
            print('Error: blastp')
            sys.exit(1)
    print('blastp done!')
    print('Done!')
    print('Output: ' + output_dir)
    print('Time: ' + str(time.time() - start_time) + 's')


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print('\nUser interrupts me! ;-) See you!')
        sys.exit(1)
