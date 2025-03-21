#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import argparse
from subprocess import Popen, PIPE
from concurrent.futures import ProcessPoolExecutor, as_completed


def shell_cmd(command) -> int:
    print(command)
    process = Popen(command, stdout=PIPE, stderr=PIPE,shell=True)

    # 等待子进程结束，并获取输出和错误信息
    return_code = process.wait()


    print(f"子进程已结束，返回码为: {return_code}")
    return return_code



def main():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--shell', type=str, help='shell command file')
    parser.add_argument('--threads', type=int, help='number of threads')
    args = parser.parse_args()
    
    shell_list = []

    shell_file = args.shell

    if not os.path.exists(shell_file):
        print('shell file not exists')
        return
    
    with open(shell_file, 'r') as f:
        for line in f.readlines():
            line = line.strip()
            if line.startswith('#'):
                continue
            shell_list.append(line)

    # 使用max_workers设置最大工作进程数为64
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        # 提交任务到进程池，并获取Future对象列表
        futures = [executor.submit(shell_cmd, shell) for shell in shell_list]
        # as_completed函数生成一个迭代器，按任务完成的顺序产生Future对象
        for future in as_completed(futures):
            # 通过Future对象获取任务结果
            result = future.result()
            print(f"任务结果: {result}")

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print('\nExiting from the program ealier!')
        sys.exit()
