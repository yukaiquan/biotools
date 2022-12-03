# psRobot_tar使用

##### 🙈: yukaiquan

##### 📧: 1962568272@qq.com

# 1.必备小知识

- psRobot_tar模块 is designed to find potential small RNA targets;
psRobot_tar 识别潜在的小RNA 的靶基因。
- psRobot_map模块 is designed to find all perfect matching locations of short sequences (less than 40bp) in longer reference sequences;
psRobot_map 在更长的参考序列上找出所有完美匹配的短序列(小于40bp)。
- psRobot_mir模块 is designed to find small RNAs with stem-loop precursors (e.g. miRNAs or shRNAs) for a batch of input sequences from high throughput sequencing data;
psRobot_mir 可为一批来自高通量的输入序列寻找具有茎环前体的小RNA(如miRNA或shRNA)。
- psRobot_deg模块 is designed to identify which small RNA targets are supported by user specified degradome data.
psRobot_deg 用于识别哪些小RNA靶标得到了用户指定的降解组数据的支持。

# 2.使用psRobot_tar识别靶基因

## 2.1准备数据
从mirBase([https://www.mirbase.org/ftp.shtml](https://www.mirbase.org/ftp.shtml))下载，mature.fa.gz 文件
## 2.2准备软件
据说依赖mfold3.5但是我没安也可以运行如果需要安装看这里

```shell
wget http://omicslab.genetics.ac.cn/psRobot/program/WebServer/psRobot_v1.2.tar.gz
tar xvzf psRobot_v1.2.tar.gz
cd  psRobot_v1.2
sudo ./configure
make
sudo make install
source  ~/.bashrc 
```
## 2.3可以从mature中挑选指定物种也可以全量预测再筛选
使用cDNA序列不用genomic序列的原因是，miRNA在细胞质和靶基因结合发挥作用。此时靶基因还有UTR区域但是已经没有内含子区了。(考虑到UTR区域的序列特点，其实用CDS序列也行)
**psRobot_tar 的参数**:

- -s input file name: smRNA sequences (fasta format)；default = smRNA
待预测的miRNA，fasta格式；默认：smRNA
- -t input file name: target sequences (fasta format)，default = target
用于搜索的cDNA序列，fasta格式；默认： target
- -o output file name，👉注意：default = smRNA-target.gTP
输出文件名，默认：smRNA-target.gTP
- -ts target penalty score, lower is better (0-5)，default = 2.5
输出结果的阈值，默认：2.5
- -fp 5 prime boundary of essential sequence (1-2)，default = 2
5'后第几位开始是必要区间(1~2)， 默认：2
- -tp 3 prime boundary of essential sequence (7-31)， default = 17
3'后第几位开始是必要区间(7~31)， 默认：17
- -gl position after which with gap/bulge permit (0-30), 0 means no gap/bulge permitted， default = 17
从第几个碱基后允许出现gap/bulge， 默认：17
- -p number of processors use，default = 1，
使用线程数， 默认：1，👉注意：根据实际情况可以改大些
- -gn number of gaps/bulges permitted (0-5)， default = 1
允许存在几个gap/bulge， 默认：1

运行吧-----
```shell
psRobot_tar -s mature.fa -t transcript.fasta -p 4 -o target_results.gTP
```
## 2.4结果查看
:::info
>dme-miR-315-5p MIMAT0000407 Drosophila melanogaster miR-315-5p Score: 2.2      LOC116195782
Query:          1 TTTTGATTGTTGCTCAGAAAGC 22
                  ||||||||*||||:||||||:*
Sbjct:       1463 AAAACTAAGAACGGGTCTTTTA 1442
>cel-miR-392-5p MIMAT0020344 Caenorhabditis elegans miR-392-5p  Score: 2.5      LOC116205351
Query:          1 AGCATTCGTGGTTGAGGATAT 21
                  *||*||||||||||||||*|*
Sbjct:        252 CCGGAAGCACCAACTCCT-TC 233
>osa-miR319a-3p.2-3p MIMAT0001028 Oryza sativa miR319a-3p.2-3p  Score: 2.0      LOC116213034
Query:          1 TTGGACTGAAGGGTGCTCCC 20
                  ||||||||||||||*|:||*
Sbjct:       2404 AACCTGACTTCCCAGGGGGA 2385
>osa-miR319b MIMAT0001029 Oryza sativa miR319b  Score: 2.0      LOC116213034
Query:          1 TTGGACTGAAGGGTGCTCCC 20
                  ||||||||||||||*|:||*
Sbjct:       2404 AACCTGACTTCCCAGGGGGA 2385
>osa-miR166k-5p MIMAT0022870 Oryza sativa miR166k-5p    Score: 2.0      LOC116205856
:::
## 2.5简化结果
```shell
grep '^>' target_results.gTP | sed -r 's#>##g' >  target_results.txt
```
:::info
dme-miR-315-5p MIMAT0000407 Drosophila melanogaster miR-315-5p  Score: 2.2      LOC116195782
cel-miR-392-5p MIMAT0020344 Caenorhabditis elegans miR-392-5p   Score: 2.5      LOC116205351
osa-miR319a-3p.2-3p MIMAT0001028 Oryza sativa miR319a-3p.2-3p   Score: 2.0      LOC116213034
osa-miR319b MIMAT0001029 Oryza sativa miR319b   Score: 2.0      LOC116213034
osa-miR166k-5p MIMAT0022870 Oryza sativa miR166k-5p     Score: 2.0      LOC116205856
osa-miR166g-5p MIMAT0022881 Oryza sativa miR166g-5p     Score: 2.5      LOC116209443
dps-miR-315 MIMAT0001261 Drosophila pseudoobscura miR-315       Score: 2.2      LOC116195782
sbi-miR319a MIMAT0001469 Sorghum bicolor miR319a        Score: 2.0      LOC116213034
ame-miR-315-5p MIMAT0001487 Apis mellifera miR-315-5p   Score: 2.2      LOC116195782
aga-miR-315 MIMAT0001521 Anopheles gambiae miR-315      Score: 2.2      LOC116195782
mghv-miR-M1-2-3p MIMAT0001565 Mouse gammaherpesvirus miR-M1-2-3p        Score: 1.5      LOC116214006
gma-miR319c MIMAT0001691 Glycine max miR319c    Score: 2.5      LOC116213034
zma-miR319a-3p MIMAT0001715 Zea mays miR319a-3p Score: 2.0      LOC116213034
zma-miR319c-3p MIMAT0001716 Zea mays miR319c-3p Score: 2.0      LOC116213034
zma-miR319b-3p MIMAT0001717 Zea mays miR319b-3p Score: 2.0      LOC116213034
zma-miR319d-3p MIMAT0001718 Zea mays miR319d-3p Score: 2.0      LOC116213034
zma-miR166k-5p MIMAT0015184 Zea mays miR166k-5p Score: 2.2      LOC116205856
ptc-miR319e MIMAT0002006 Populus trichocarpa miR319e    Score: 2.5      LOC116213034
ptc-miR319f MIMAT0002007 Populus trichocarpa miR319f    Score: 2.5      LOC116213034
ptc-miR319g MIMAT0002008 Populus trichocarpa miR319g    Score: 2.5      LOC116213034
ptc-miR319h MIMAT0002009 Populus trichocarpa miR319h    Score: 2.5      LOC116213034
hsa-miR-574-5p MIMAT0004795 Homo sapiens miR-574-5p     Score: 1.5      LOC116195782
mmu-miR-694 MIMAT0003474 Mus musculus miR-694   Score: 2.5      LOC116203787
xtr-miR-29c-5p MIMAT0003676 Xenopus tropicalis miR-29c-5p       Score: 2.5      LOC116187206
:::

# 引用：

[PsRobot的使用流程 | 小麦篇](https://www.jianshu.com/p/735a89879f39)

[miRBase](https://www.mirbase.org/search.shtml)
