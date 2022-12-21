#!/usr/bin/env python3
import sys
from subprocess import Popen, PIPE


def run_fastp(fastq1, fastq2, outdir, threads):
    cmd = f'fastp -i {fastq1} -I {fastq2} -o {outdir}/fastp_out1.fq.gz -O {outdir}/fastp_out2.fq.gz -h {outdir}/fastp.html -j {outdir}/fastp.json -w {threads}'
    p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    stdout, stderr = p.communicate()
    print(stdout)
    print(stderr)
    return stdout, stderr


def main():
    fastq1 = '/home/yukaiquan/Ex-seq/data/ERR030061_1.fastq.gz'
    fastq2 = '/home/yukaiquan/Ex-seq/data/ERR030061_2.fastq.gz'
    outdir = '/home/yukaiquan/Ex-seq/data'
    threads = 16
    run_fastp(fastq1, fastq2, outdir, threads)


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print('KeyboardInterrupt')
        sys.exit(0)
