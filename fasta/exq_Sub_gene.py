'''
Author: yukaiquan
Date: 2025-05-21 21:31:15
Email: 1962568272@qq.com
LastEditor: yukaiquan
'''
import re
import os
import sys
import argparse
import gzip
from tqdm import tqdm


def main():
    args = check_options(get_options())
    fasta_file = args.input
    out_file = args.output
    gff_file = args.gff
    gene_2_chr = {}
    chr_name_list = []
    # with open(gff_file, 'r') as f:
    with gzip.open(gff_file, 'rt') as f:
        for line in tqdm(f, desc="Reading gff file..."):
            l = line.strip()
            if l.startswith("#"):
                continue
            lines = l.split("\t")
            chr_name = lines[0]
            # 自动大小写chrun cpDNA ChrUnknown ChrSy 不区分大小写
            if chr_name.lower() == "cpdna" or chr_name.lower() == "chrunknown" or chr_name.lower() == "chrsy":
                continue
            type_name = lines[2]
            if type_name.lower() != "gene":
                continue
            gene_name = lines[8].split(";")[0].split("=")[1].strip()
            if gene_name not in gene_2_chr:
                gene_2_chr[gene_name] = chr_name
            else:
                print("Warning: %s is in %s and %s" %
                      (gene_name, gene_2_chr[gene_name], chr_name))
                gene_2_chr[gene_name] = chr_name
            if chr_name not in chr_name_list:
                chr_name_list.append(chr_name)
    # 拆分亚基因组例如：chr1A chr1B chr1C chr1D 拆分为A B C D
    # chr_name_list = [re.sub(r"^.*?(\D+)$", r"\1", chr_name)
    #                  for chr_name in chr_name_list]
    chr_name_2_sub = {}
    for chr_name in chr_name_list:
        if re.search(r"^.*?(\D+)$", chr_name):
            sub_name = re.search(r"^.*?(\D+)$", chr_name).group(1)
            if chr_name not in chr_name_2_sub:
                chr_name_2_sub[chr_name] = sub_name
            else:
                print("Warning: %s is in %s and %s" %
                      (chr_name, chr_name_2_sub[chr_name], sub_name))
                chr_name_2_sub[chr_name] = sub_name
    # 处理序列中的特殊字符串
    p = re.compile(r'[-,$()#+&*.]')
    # 原始序列哈希表
    sequence = {}
    with open(fasta_file, 'r') as f:
        for line in tqdm(f, desc="Reading fasta file..."):
            l = line.strip()
            if l.startswith(">"):
                name = l[1:]
                name = re.sub(r'[ \t].*$', "", name)
                sequence[name] = ''
            else:
                l = re.sub(p, '', l)
                # 序列转大写
                sequence[name] += l.upper()
    # 根据不同的sub_name进行创建不同的文件
    for chr_name in chr_name_2_sub:
        sub_name = chr_name_2_sub[chr_name]
        file_name = out_file + "_" + sub_name + ".fasta"
        # with open(file_name, 'w') as f:
        # 如果文件不存在则创建文件 如果存在则追加内容
        if not os.path.exists(file_name):
            with open(file_name, 'w') as f:
                for gene_name in gene_2_chr:
                    if gene_2_chr[gene_name] == chr_name:
                        if gene_name not in sequence:
                            print("Warning: %s is not in %s" %
                                  (gene_name, fasta_file))
                            continue
                        seq = sequence[gene_name]
                        f.write('>' + gene_name + '\n')
                        # 更改每行碱基的长度为60
                        result1 = [seq[i:i + 60]
                                   for i in range(0, len(seq), 60)]
                        result = '\n'.join(result1)
                        f.write(result + '\n')
                        print("Writing %s to %s" % (gene_name, file_name))
        else:
            with open(file_name, 'a') as f:
                for gene_name in gene_2_chr:
                    if gene_2_chr[gene_name] == chr_name:
                        if gene_name not in sequence:
                            print("Warning: %s is not in %s" %
                                  (gene_name, fasta_file))
                            continue
                        seq = sequence[gene_name]
                        f.write('>' + gene_name + '\n')
                        # 更改每行碱基的长度为60
                        result1 = [seq[i:i + 60]
                                   for i in range(0, len(seq), 60)]
                        result = '\n'.join(result1)
                        f.write(result + '\n')
                        print("Writing %s to %s" % (gene_name, file_name))
    print("All done!")
    


def get_options():
    print('##############################################################')
    print('# yukaiquan 1962568272@qq.com')
    print('##############################################################')
    parser = argparse.ArgumentParser(
        description="Find Sub transcript by pep and gff. author: yukaiquan email: 1962568272@qq.com",
        prog='Ex_Sub_pep',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Example:\n"
        "python exq_Sub_gene.py -i fasta.pep -g gene.gff3 -o sfs  \\ \n"
    )

    parser.add_argument('-V', '--version', action='version',
                        version='%(prog)s 1.0')

    parser.add_argument(
        '-i',
        '--input',
        dest='input',
        help='longest pep fasta file',
        required=True)

    parser.add_argument(
        '-o',
        '--output',
        dest='output',
        help='The output files prefix',
        required=True)

    parser.add_argument(
        '-g',
        '--gff',
        dest='gff',
        help='gene gff file',
        required=True)

    return parser

def check_options(parser):

    args = parser.parse_args()

    if not os.path.exists(args.input):

        print("Can not locate input fasta protain, please input full path of protain\n")

        parser.print_help()

        sys.exit(1)

    if not os.path.exists(args.gff):
            print("Can not locate gff file, beacause it exists!\n")
            parser.print_help()
            sys.exit(1)

    return args


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(e)
        sys.exit(1)