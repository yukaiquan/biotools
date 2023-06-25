import sys

# Read in the blast output file
blastout = sys.argv[1]
out_file = sys.argv[2]

chr_name = ['chr1A', 'chr2A', 'chr3A', 'chr4A', 'chr5A', 'chr6A', 'chr7A', 'chr1C', 'chr2C', 'chr3C',
            'chr4C', 'chr5C', 'chr6C', 'chr7C', 'chr1D', 'chr2D', 'chr3D', 'chr4D', 'chr5D', 'chr6D', 'chr7D', 'chrUn']

# 统计列表中元素出现的次数


def count_list(list):
    dict = {}
    for i in list:
        if i not in dict:
            dict[i] = 1
        else:
            dict[i] += 1
    return dict


id_dict = {}
with open(blastout, 'r') as r:
    for line in r:
        line = line.strip().split('\t')
        if line[0] not in id_dict:
            id_dict[line[0]] = [line[1]]
        else:
            id_dict[line[0]].append(line[1])

with open(out_file, 'w') as w:
    for key in id_dict:
        id_dict[key] = count_list(id_dict[key])
        w.write(key+'\t')
        for chr in chr_name:
            if chr not in id_dict[key]:
                w.write('0'+'\t')
            else:
                w.write(str(id_dict[key][chr])+'\t')
        w.write('\n')
