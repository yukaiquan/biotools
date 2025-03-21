import sys
import os

def TidyMirandaResult(inputfile, outfile)->dict:
    out_dict = {}
    infpath = inputfile
    # 获取inputfile的路径
    outfpath = os.path.join(os.path.dirname(infpath), outfile)
    output_open = open(outfpath, 'w')
    with open(infpath,'r') as f:
        for row in f:
            row = row.split('\t')
            if ',' in row[0]:
                mi_name = row[0][2:].split(',')[0]
            else:
                mi_name = row[0][2:]
            if ']' in mi_name:
                mi_name = mi_name.replace('[', '')
            if ']' in mi_name:
                mi_name = mi_name.replace(']', '')
            gene_name = row[1]
            if mi_name not in out_dict.keys():
                out_dict[mi_name] = [gene_name]
            else:
                out_dict[mi_name].append(gene_name)
            # 写入文件
            line = mi_name + 't' + gene_name + '\n'
            output_open.writelines(line)
    output_open.close()
    return out_dict

def TidyTargetscanResult(path, inputfile, outfile):
    infpath = r'{}\{}'.format(path, inputfile)
    outfpath = r'{}\{}'.format(path, outfile)
    i = 0
    with open(infpath) as f:
        for row in f:
            i += 1
            if i == 1:
                continue
            row = row.split('\t')
            line = row[1] + ':' + row[0] + '\n'
            with open(outfpath, 'a') as r:
                r.writelines(line)

def TidyRNA22Result(path, inputfile, outfile):
    infpath = r'{}\{}'.format(path, inputfile)
    outfpath = r'{}\{}'.format(path, outfile)
    with open(infpath) as f:
        for row in f:
            row = row.split('\t')
            line = row[0] + ':' + row[1] + '\n'
            with open(outfpath, 'a') as r:
                r.writelines(line)

def targetfinder(inputfile, outfile)->dict:
    out_dict = {}
    infpath = inputfile
    outfpath = os.path.join(os.path.dirname(infpath), outfile)
    out_file_open = open(outfpath, 'w')
    with open(infpath, 'r') as f:
        for row in f:
            if row.startswith('query='):
                lines = row.strip().split(',')
                mi_name = lines[0].split('=')[1].strip()
                if ']' in mi_name:
                    mi_name = mi_name.replace('[', '')
                if ']' in mi_name:
                    mi_name = mi_name.replace(']', '')
                # 查找target在list的位置
                target_index = [i for i, x in enumerate(lines) if 'target=' in x]
                if len(target_index) == 0:
                    continue
                else:
                    target_index = target_index[0]
                    gene_name = lines[target_index].split('=')[1].strip()

                if mi_name not in out_dict.keys():
                    out_dict[mi_name] = [gene_name]
                else:
                    out_dict[mi_name].append(gene_name)
                out_file_open.write(mi_name + ':' + gene_name + '\n')
    out_file_open.close()    
    
    return out_dict

def targetpsRoutine(inputfile, outfile)->dict:
    out_dict = {}
    infpath = inputfile
    outfpath = os.path.join(os.path.dirname(infpath), outfile)
    out_file_open = open(outfpath, 'w')
    with open(infpath, 'r') as f:
        for row in f:
            if row.startswith('>'):
                lines = row.strip().split('\t')
                mi_name = lines[0].replace('>', '')
                if ']' in mi_name:
                    mi_name = mi_name.replace('[', '')
                if ']' in mi_name:
                    mi_name = mi_name.replace(']', '')
                gene_name = lines[2]
                if mi_name not in out_dict.keys():
                    out_dict[mi_name] = [gene_name]
                else:
                    out_dict[mi_name].append(gene_name)
                out_file_open.write(mi_name + ':' + gene_name + '\n')
    out_file_open.close()
    return out_dict


# 查找2个dict的交集
def find_intersection(dict1, dict2)->dict:
    intersection = {}
    for key in dict1.keys():
        value_list1 = dict1[key]
        try:
            value_list2 = dict2[key]
        except:
            print(key)
            continue
        for value in value_list1:
            if value in value_list2:
                if key not in intersection.keys():
                    intersection[key] = [value]
                else:
                    intersection[key].append(value)

    return intersection


def main():
    miranda = sys.argv[1]
    finder = sys.argv[2]
    targetpsRobot = sys.argv[3]
    merge = sys.argv[4]
    deg_list = sys.argv[5]
    deg_out = sys.argv[6]
    # path = r'D:\用户\桌面\练习\结果'
    miranda_dict = TidyMirandaResult(miranda, 'miranda_TargetResult.txt')
    finder_dict = targetfinder(finder, 'targetfinder_TargetResult.txt')
    intersection = find_intersection(miranda_dict, finder_dict)
    targetpsRobot_dict = targetpsRoutine(targetpsRobot, 'targetpsRobots_TargetResult.txt')
    intersections = find_intersection(intersection, targetpsRobot_dict)
    with open(merge, 'w') as f:
        for key in intersections.keys():
            value_list = intersections[key]
            for value in value_list:
                line = key + ':' + value + '\n'
                f.writelines(line)
    
    # 处理deg_list
    deg = []
    with open(deg_list, 'r') as f:
        for row in f:
            deg.append(row.strip())
    deg_out_open = open(deg_out, 'w')
    for key in intersection.keys():
        for value in deg:
            if value in key:
                value_list = intersection[key]
                for value in value_list:
                    line = key + ':' + value + '\n'
                    deg_out_open.write(line)
    deg_out_open.close()


    # TidyRNA22Result(path, 'RNA22_result.txt', 'RNA22_TidyResult.txt')
    # TidyTargetscanResult(path, 'targetscan_result.txt', 'targetscan_TidyResult.txt')
    # TidyTargetscanResult(path, 'pita_results.txt', 'pita_TidyResult.txt') # pita结果处理和targetscan是一样的



if __name__ == '__main__':
    main()