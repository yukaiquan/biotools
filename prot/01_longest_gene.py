'''
Author: yukaiquan
Date: 2022-07-26 16:31:15
Email: 1962568272@qq.com
LastEditor: yukaiquan
Description: 01_longest_gene.py
'''
import re
import sys

def main():
    args = check_options(get_options())
    fasta_file = args.input
    print("Reading fasta file...")
    out_file = args.output
    # 原始序列哈希表
    sequence = {}
    with open(fasta_file, 'r') as f:
        for line in f:
            l = line.strip()
            if l.startswith(">"):
                name = l[1:]
                name = re.sub(" .*$", "", name)
                sequence[name] = ''
            else:
                sequence[name] += l
    transcrip_length = {}
    # 获取每个转录本长度 并存入哈希表
    for name, seq in sequence.items():
        transcrip_length[name] = len(seq)
    # 最长转录本哈希表
    longest_transcrip = {}
    # 定义转录本名称的正则表达式 转录本名称的格式可以为 AVESA.00010b.r2.1AG0069130.RA AVESA.00010b.r2.1AG0069130.2 AVESA.00010b.r2.1AG0069130-RA AVESA.00010b.r2.1AG0069130-2
    pattarn = re.compile(r"^.*\d+(?=\.\w+$|\.\d+$|-\w+$|-\d+$)")
    for name, length in transcrip_length.items():
        gene_name = re.search(pattarn, name).group()
        if gene_name not in longest_transcrip:
            longest_transcrip[gene_name] = name
        else:
            if length > transcrip_length[name]:
                longest_transcrip[gene_name] = name    # 更新最长转录本
    # 提取序列生成文件
    print("write fasta file...")
    with open(out_file, 'w') as f:
        for gene_name, trans_name in longest_transcrip.items():
            f.write('>' + gene_name + '\n')
            # 更改每行碱基的长度为60
            result1 = [sequence[trans_name][i:i + 60]
                       for i in range(0, len(sequence[trans_name]), 60)]
            result = '\n'.join(result1)
            f.write(result + '\n')

def get_options():

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

    return args


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(e)
        sys.exit(1)