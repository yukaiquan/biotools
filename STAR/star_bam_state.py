#!/use/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import click


# 列出目录下后缀为Log.final.out的文件
def listdir(path: str, suffix='Log.final.out') -> list:
    list_file = []
    # 只需要当前目录
    for file in os.listdir(path):
        if os.path.isfile(os.path.join(path, file)) and file.endswith(suffix):
            list_file.append(os.path.join(path, file))
    return list_file


def analysis_log(bam_out_list: list, out_file: str) -> bool:
    out_file_f = open(out_file, 'w')
    out_file_f.write('file\tinput reads\tuni mapped reads\tuni mapped percent\tmulti mapped reads\tmulti mapped percent\ttoo many mapped reads\ttoo many mapped percent\ttoo many mis reads\ttoo many mis percent\ttoo short mis reads\ttoo short mis percent\tmis other\tmis other percent\n')

    # print(bam_out_list)
    for log_file in bam_out_list:
        with open(log_file, 'r') as f:
            for line in f:
                # print(line)
                if 'Number of input reads' in line:
                    input_reads = line.split('|')[1].strip()
                elif 'Uniquely mapped reads number' in line:
                    uni_mapped_reads = line.split('|')[1].strip()
                elif 'Uniquely mapped reads %' in line:
                    uni_mapped_reads_percent = line.split('|')[1].strip()
                # Multi-mapped reads
                elif 'Number of reads mapped to multiple loci' in line:
                    multi_mapped_reads = line.split('|')[1].strip()
                elif '% of reads mapped to multiple loci' in line:
                    multi_mapped_reads_percent = line.split('|')[1].strip()
                elif 'Number of reads mapped to too many loci' in line:
                    too_mapped_reads = line.split('|')[1].strip()
                elif '% of reads mapped to too many loci' in line:
                    too_mapped_reads_percent = line.split('|')[1].strip()
                # Unmapped reads
                elif 'Number of reads unmapped: too many mismatches' in line:
                    too_many_mismatches = line.split('|')[1].strip()
                elif '% of reads unmapped: too many mismatches' in line:
                    too_many_mismatches_percent = line.split('|')[1].strip()
                elif 'Number of reads unmapped: too short' in line:
                    too_short = line.split('|')[1].strip()
                elif '% of reads unmapped: too short' in line:
                    too_short_percent = line.split('|')[1].strip()
                elif 'Number of reads unmapped: other' in line:
                    other = line.split('|')[1].strip()
                elif '% of reads unmapped: other' in line:
                    other_percent = line.split('|')[1].strip()
                else:
                    continue
            out_file_f.write(log_file + '\t' + input_reads + '\t' + uni_mapped_reads + '\t' + uni_mapped_reads_percent + '\t' + multi_mapped_reads + '\t' + multi_mapped_reads_percent + '\t' + too_mapped_reads +
                             '\t' + too_mapped_reads_percent + '\t' + too_many_mismatches + '\t' + too_many_mismatches_percent + '\t' + too_short + '\t' + too_short_percent + '\t' + other + '\t' + other_percent + '\n')
    out_file_f.close()
    return True

# 命令行提示
@click.command(name='star_bam_state.py', help='analysis star bam state')
@click.version_option(version='1.0.0', prog_name='star_bam_state.py', message='%(prog)s version %(version)s')
@click.option('--bam_path', '-b', help='star output bam file path')
@click.option('--out_file', '-o', help='out file')
# 自动调帮助文档
@click.help_option('-h', '--help', help='show this help message and exit')
def main(bam_path, out_file):
    bam_out_list = listdir(bam_path, 'Log.final.out')
    if analysis_log(bam_out_list, out_file):
        print('analysis success')
    else:
        print('analysis failed')


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print('Interrupted')
        sys.exit(0)
