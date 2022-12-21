##### 🙈: yukaiquan
##### 📧: 1962568272@qq.com
# 1 STAR的比对结果
储存在后缀Log.final.out的文件中
>                                  Started job on |	Dec 09 23:50:24
>                              Started mapping on |	Dec 10 00:08:56
>                                     Finished on |	Dec 10 00:16:18
>        Mapping speed, Million of reads per hour |	208.33
> 
>                           Number of input reads |	25578796
>                       Average input read length |	298
>                                     UNIQUE READS:
>                    Uniquely mapped reads number |	21730084
>                         Uniquely mapped reads % |	84.95%
>                           Average mapped length |	293.45
>                        Number of splices: Total |	19257323
>             Number of splices: Annotated (sjdb) |	19145871
>                        Number of splices: GT/AG |	18746377
>                        Number of splices: GC/AG |	238662
>                        Number of splices: AT/AC |	8581
>                Number of splices: Non-canonical |	263703
>                       Mismatch rate per base, % |	0.20%
>                          Deletion rate per base |	0.00%
>                         Deletion average length |	1.00
>                         Insertion rate per base |	0.01%
>                        Insertion average length |	2.65
>                              MULTI-MAPPING READS:
>         Number of reads mapped to multiple loci |	2637165
>              % of reads mapped to multiple loci |	10.31%
>         Number of reads mapped to too many loci |	40102
>              % of reads mapped to too many loci |	0.16%
>                                   UNMAPPED READS:
>   Number of reads unmapped: too many mismatches |	0
>        % of reads unmapped: too many mismatches |	0.00%
>             Number of reads unmapped: too short |	1152759
>                  % of reads unmapped: too short |	4.51%
>                 Number of reads unmapped: other |	18686
>                      % of reads unmapped: other |	0.07%
>                                   CHIMERIC READS:
>                        Number of chimeric reads |	0
>                             % of chimeric reads |	0.00%

# 2 有很多样怎么一次性统计
那就随便写一个脚本吧
```python
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

```
脚本github仓库
# 3 结果
![image.png](https://cdn.nlark.com/yuque/0/2022/png/32371416/1671633901613-017fc330-0573-4475-8b1c-094dd7bab3ca.png#averageHue=%23efefef&clientId=u6923d932-ae1f-4&crop=0&crop=0&crop=1&crop=1&from=paste&height=710&id=u716c4674&margin=%5Bobject%20Object%5D&name=image.png&originHeight=710&originWidth=1814&originalType=binary&ratio=1&rotation=0&showTitle=false&size=227257&status=done&style=none&taskId=u6fffea8e-f1b7-4ef3-8aa0-dd3fc90faf5&title=&width=1814)

