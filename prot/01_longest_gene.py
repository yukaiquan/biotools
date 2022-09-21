'''
Author: yukaiquan
Date: 2022-07-26 16:31:15
Email: 1962568272@qq.com
LastEditor: yukaiquan
Description: python /mnt/e/software/Ex-seq/biotools/prot/01_longest_gene.py -i /mnt/d/databasezip/genome/ACD/Sang/gene_annotation_1/gene_annotation/ACD_sativa_sang_cds.fasta -o test.fa
'''
import re
import os
import sys
import argparse
from tqdm import tqdm


def main():
    args = check_options(get_options())
    fasta_file = args.input
    # fasta_file = r"D:\databasezip\genome\ACD\Sang\gene_annotation_1\gene_annotation\ACD_sativa_sang_cds.fasta"
    out_file = args.output
    # out_file = "E:\\oatdatabase\\11_genefamily\\01_longest_pep\\test.fasta"
    # 处理序列中的特殊字符串
    p = re.compile(r'[-,$()#+&*]')
    # 原始序列哈希表
    sequence = {}
    with open(fasta_file, 'r') as f:
        for line in tqdm(f, desc="Reading fasta file..."):
            l = line.strip()
            if l.startswith(">"):
                name = l[1:]
                name = re.sub(" .*$", "", name)
                sequence[name] = ''
            else:
                l = re.sub(p, '', l)
                # 序列转大写
                sequence[name] += l.upper()
    transcrip_length = {}
    # 获取每个转录本长度 并存入哈希表
    for name, seq in sequence.items():
        transcrip_length[name] = len(seq)
    # 最长转录本哈希表
    longest_transcrip = {}
    # 定义转录本名称的正则表达式 转录本名称的格式可以为 AVESA.00010b.r2.1AG0069130.RA AVESA.00010b.r2.1AG0069130.2 AVESA.00010b.r2.1AG0069130-RA AVESA.00010b.r2.1AG0069130-2
    # 20220921 add AVESA.00010b.r2.1AG0069130LC.1
    pattarn = re.compile(r"^.*\w+(?=\.\w+$|\.\d+$|-\w+$|-\d+$)")
    for name, length in tqdm(transcrip_length.items(), desc="Finding longest transcript..."):
        gene_name = re.search(pattarn, name).group()
        if gene_name not in longest_transcrip:
            longest_transcrip[gene_name] = name
        else:
            if length > transcrip_length[name]:
                print("Warning: %s is longer than %s" %
                      (name, longest_transcrip[gene_name]))
                longest_transcrip[gene_name] = name    # 更新最长转录本
    if args.gene_2_transcript:
        with open(args.gene_2_transcript, 'w') as f:
            for gene, transcript in longest_transcrip.items():
                f.write(gene + "\t" + transcript + "\n")
    # 提取序列生成文件
    with open(out_file, 'w') as f:
        for gene_name, trans_name in tqdm(longest_transcrip.items(), desc="Writing longest transcript to file..."):
            f.write('>' + gene_name + '\n')
            # 更改每行碱基的长度为60
            result1 = [sequence[trans_name][i:i + 60]
                       for i in range(0, len(sequence[trans_name]), 60)]
            result = '\n'.join(result1)
            f.write(result + '\n')


def get_options():
    print('##############################################################')
    print('# yukaiquan 1962568272@qq.com')
    print('##############################################################')
    parser = argparse.ArgumentParser(
        description="Find longest transcript inded gene pep (sequence name splite with '.'). author: yukaiquan email: 1962568272@qq.com",
        prog='Ex_Longest_pep',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Example:\n"
        "python 03_Ex_longest_tr.py -i fasta.pep -o long_seq.fasta  \\ \n"
    )

    parser.add_argument('-V', '--version', action='version',
                        version='%(prog)s 1.0')

    parser.add_argument(
        '-i',
        '--input',
        dest='input',
        help='transcrit pep fasta file',
        required=True)

    parser.add_argument(
        '-o',
        '--output',
        dest='output',
        help='The output files for saving results',
        required=True)

    parser.add_argument(
        '-g',
        '--gene_2_transcript',
        dest='gene_2_transcript',
        help='longest transcript from where file(default: None)',
        required=False)

    return parser


def check_options(parser):

    args = parser.parse_args()

    if not os.path.exists(args.input):

        print("Can not locate input fasta protain, please input full path of protain\n")

        parser.print_help()

        sys.exit(1)

    # End check bwa

    if not args.output:

        print("Can not locate output file name, please input output file name.\n")

        parser.print_help()

        sys.exit(1)
    if args.gene_2_transcript:
        if os.path.exists(args.gene_2_transcript):
            print("Can not locate gene_2_transcript file, beacause it exists!\n")
            parser.print_help()
            sys.exit(1)

    return args


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(e)
        sys.exit(1)
