import numpy as np
import sys
import tqdm as tqdm

input_file = sys.argv[1]
# chr1A  1000000 0.3
# chr1A  2000000 0.5
# chr1A  2300000 0.3
# 2400000 0.3
# 2600000 0.3
# 3000000 0.3
# 4000000 0.2
# 1000000 0.3
# 2000000 0.5
# 2300000 0.3
# 2400000 0.3
# 2600000 0.3
# 指定染色体长度bp
input_chr_len_file = sys.argv[2]
# 指定滑窗大小1Mb
window_size = int(sys.argv[3])
output_file = sys.argv[4]
# output_png = sys.argv[5]

# pos_dict = {1000000: 0.3, 2000000: 0.5, 2300000: 0.3, 2400000: 0.3, 2600000: 0.3,
#             3000000: 0.3, 4000000: 0.2}


def main():
    with open(input_chr_len_file, 'r') as f:
        lines = f.readlines()
        chr_len_dict = {}
        for line in lines:
            line = line.strip().split('\t')
            chr_len_dict[line[0]] = int(line[1])
    output = open(output_file, 'w')
    # 高效迭代字典
    for chr, value in tqdm.tqdm(chr_len_dict.items(), desc='开始处理染色体'):
        print('#'*25 + chr + '#'*25)
        pos_dict = {}
        with open(input_file, 'r') as f:
            lines = f.readlines()
            for line in tqdm.tqdm(lines, desc='读取数据'):
                line = line.strip().split('\t')
                pos_dict[int(line[1])] = float(line[2])
        chr_lenth = generate_interval_list(value, window_size)
        chr_dict = generate_interval_dict(chr_lenth)

        for key, value in tqdm.tqdm(chr_dict.items(), desc='分配区间'):
            key = key.split('-')
            for pos, cov in pos_dict.items():
                # 起始位置在区间内 终止位置不在区间内 但是如果终止位置等于染色体长度的话，也要加上
                if pos >= int(key[0]) and pos < int(key[1]):
                    value.append(cov)
        for key, value in tqdm.tqdm(chr_dict.items(), desc='计算平均值'):
            # 如果区间内没有数据，那么就赋值为0
            if len(value) == 0:
                chr_dict[key] = 0
                continue
            value = np.array(value, dtype=float)
            value = np.mean(value)
            chr_dict[key] = value
        for key, value in tqdm.tqdm(chr_dict.items(), desc='生成数据'):
            key = key.split('-')
            pos = int(key[0])+(int(key[1]) - int(key[0]))/2
            output.write(chr + '\t' + str(int(pos)) + '\t' + str(value) + '\n')
    output.close()

    # {'0-1000000': 0,
    #  '1000000-2000000': 0.3,
    #  '2000000-3000000': 0.35000000000000003,
    #  '3000000-4000000': 0.3,
    #  '4000000-5000000': 0.2}
    # 根据chr_dict生成x轴和y轴的数据

    # x = []
    # # y = []
    # output = open(output_file, 'w')
    # for key, value in tqdm.tqdm(chr_dict.items(), desc='生成数据'):
    #     key = key.split('-')
    #     pos = int(key[0])+(int(key[1]) - int(key[0]))/2
    #     output.write(str(int(pos)) + '\t' + str(value) + '\n')
    # output.close()


# 利用numpy生成区间列表


def generate_interval_list(genome_chr_len, window_size) -> list:
    # interval_list = np.arange(0, genome_chr_len + 1, window_size)
    # 最后一个区间的长度可能不足1Mb，所以需要单独处理
    interval_list = np.arange(0, genome_chr_len, window_size)
    interval_list = np.append(interval_list, genome_chr_len)
    interval_list = interval_list.astype(str)
    return interval_list

# 根据区间列表生成区间字典


def generate_interval_dict(interval_list) -> dict:
    interval_dict = {}
    for i in range(len(interval_list) - 1):
        interval_dict[interval_list[i] + '-' + interval_list[i + 1]] = []
    return interval_dict


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me! ;-) Bye!")
