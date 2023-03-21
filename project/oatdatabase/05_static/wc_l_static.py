# -*- coding:utf-8 -*-
import os
import sys

ROOT_PATH = sys.argv[1]
COUNT_TYPE = ['tsx', 'ts', 'js', 'scss', 'html', 'css']
if __name__ == '__main__':
    lines = 0
    for filepath, dirnames, filenames in os.walk(ROOT_PATH):
        for filename in filenames:
            path = os.path.join(filepath, filename)
            # 跳过.开头的文件夹
            if '/.' in path:
                continue
            # 去除没有后缀的文件
            elif len(filename.split(".")) == 1:
                continue
            # 跳过node_modules文件夹
            elif 'node_modules' in path:
                continue
            else:
                type = filename.split(".")[1]
                if(type in COUNT_TYPE):
                    count = len(open(path, encoding='UTF-8').readlines())
                    print(path, ",lines:", count)
                    lines += count
    print("total count :", lines)
