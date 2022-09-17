##### 🙈: yukaiquan
##### 📧: 1962568272@qq.com
# 整理染色体编号
## 更改零碎的片段
存在76条染色体序列，实际07之后只是一些小片段，为了方便数据库存储，更改07之后的序列以100N相连组成chrUn序列
> Chr01   583880110       7       60      61
> Chr04   509801500       593611460       509801500       509801501
> Chr06   462827500       1103412968      462827500       462827501
> Chr07   452504500       1566240476      452504500       452504501
> Chr02   523817500       2018744984      523817500       523817501
> Chr05   475748926       2542562492      475748926       475748927
> Chr03   518309500       3018311426      518309500       518309501
> Chr08   42000   3536620934      42000   42001
> Chr09   37520   3536662942      37520   37521
> Chr10   133592  3536700470      133592  133593
> Chr11   32385   3536834070      32385   32386
> Chr12   28354   3536866463      28354   28355
> Chr13   28246   3536894825      28246   28247
> Chr14   80000   3536923079      80000   80001
> Chr15   44000   3537003087      44000   44001
> Chr16   103000  3537047095      103000  103001
> Chr17   344000  3537150103      344000  344001
> Chr18   130500  3537494111      130500  130501
> Chr19   62000   3537624619      62000   62001
> Chr20   43000   3537686627      43000   43001
> Chr21   64448   3537729635      64448   64449
> Chr22   375000  3537794091      375000  375001
> Chr23   646000  3538169099      646000  646001
> Chr24   53000   3538815107      53000   53001
> .......

使用脚本合并该部分序列
`python3 /mnt/e/python3/fastaN60/fastaN60.py -l contig.list`
contig.list文件按照该脚本的顺序进行合并：
> Chr08.fasta
> Chr09.fasta
> Chr10.fasta
> Chr11.fasta
> Chr12.fasta
> ....

脚本链接：
## 更改01-07的染色体编号
| Chr01 | chr1A | 583880110 |
| --- | --- | --- |
| Chr02 | chr2A | 523817500 |
| Chr03 | chr3A | 518309500 |
| Chr04 | chr4A | 509801500 |
| Chr05 | chr5A | 475748926 |
| Chr06 | chr6A | 462827500 |
| Chr07 | chr7A | 452504500 |

# 更改cds和pep的ID
因为每个基因只有一个转录本所以只需要在原来的cds和pep文件后方加.1就好同时去掉protain后面的*

```bash
awk '{if (/^>/) {gsub(/$/, ".1"); print} else if (/^.*$/){ print }}' A.strigosa.chr.HCC.cds > A_strigosa_nc_cds.fasta
awk '{if (/^>/) {gsub(/$/, ".1"); print} else if (/^.*$/){ print }}' ../A.strigosa.chr.HCC.pep | sed 's#\*##g' > A_strigosa_nc_pep.fasta
```
## 更改GFF文件
开始我还以为GFF文件没什么大问题直到我看到这里
> chr1A   EVM     gene    369603  371087  .       -       .       ID=AS01G000010;Name=AS01G000010
> chr1A   EVM     mRNA    369603  371087  .       -       .    ID=AS01G000010;Parent=AS01G000010;Name=AS01G000010
> chr1A   EVM     exon    369603  371087  .       -       .       ID=AS01G000010.exon1;Parent=AS01G000010
> chr1A   EVM     CDS     369603  371087  .       -       0       ID=cds.AS01G000010;Parent=AS01G000010

**mRNA**的parent ID应为gene的ID
**exon cds**的parent ID应为mRNA的ID
同时还存在以下问题cds编号并没有像exon的编号在变化
> chr1A   EVM     gene    573979  576937  .       +       .       ID=AS01G000030;Name=AS01G000030
> chr1A   EVM     mRNA    573979  576937  .       +       .       ID=AS01G000030;Parent=AS01G000030;Name=AS01G000030
> chr1A   EVM     exon    573979  574093  .       +       .       ID=AS01G000030.exon1;Parent=AS01G000030
> chr1A   EVM     CDS     573979  574093  .       +       0       ID=cds.AS01G000030;Parent=AS01G000030
> chr1A   EVM     exon    574176  574270  .       +       .       ID=AS01G000030.exon2;Parent=AS01G000030
> chr1A   EVM     CDS     574176  574270  .       +       2       ID=cds.AS01G000030;Parent=AS01G000030
> chr1A   EVM     exon    574381  574472  .       +       .       ID=AS01G000030.exon3;Parent=AS01G000030
> chr1A   EVM     CDS     574381  574472  .       +       0       ID=cds.AS01G000030;Parent=AS01G000030
> chr1A   EVM     exon    574969  575062  .       +       .       ID=AS01G000030.exon4;Parent=AS01G000030
> chr1A   EVM     CDS     574969  575062  .       +       1       ID=cds.AS01G000030;Parent=AS01G000030
> chr1A   EVM     exon    575143  575215  .       +       .       ID=AS01G000030.exon5;Parent=AS01G000030
> chr1A   EVM     CDS     575143  575215  .       +       0       ID=cds.AS01G000030;Parent=AS01G000030
> chr1A   EVM     exon    575304  575436  .       +       .       ID=AS01G000030.exon6;Parent=AS01G000030
> chr1A   EVM     CDS     575304  575436  .       +       2       ID=cds.AS01G000030;Parent=AS01G000030
> chr1A   EVM     exon    575572  575746  .       +       .       ID=AS01G000030.exon7;Parent=AS01G000030
> chr1A   EVM     CDS     575572  575746  .       +       1       ID=cds.AS01G000030;Parent=AS01G000030
> chr1A   EVM     exon    575907  576047  .       +       .       ID=AS01G000030.exon8;Parent=AS01G000030
> chr1A   EVM     CDS     575907  576047  .       +       0       ID=cds.AS01G000030;Parent=AS01G000030
> chr1A   EVM     exon    576117  576212  .       +       .       ID=AS01G000030.exon9;Parent=AS01G000030
> chr1A   EVM     CDS     576117  576212  .       +       0       ID=cds.AS01G000030;Parent=AS01G000030
> chr1A   EVM     exon    576318  576459  .       +       .       ID=AS01G000030.exon10;Parent=AS01G000030
> chr1A   EVM     CDS     576318  576459  .       +       0       ID=cds.AS01G000030;Parent=AS01G000030
> chr1A   EVM     exon    576552  576679  .       +       .       ID=AS01G000030.exon11;Parent=AS01G000030
> chr1A   EVM     CDS     576552  576679  .       +       2       ID=cds.AS01G000030;Parent=AS01G000030
> chr1A   EVM     exon    576761  576937  .       +       .       ID=AS01G000030.exon12;Parent=AS01G000030
> chr1A   EVM     CDS     576761  576937  .       +       0       ID=cds.AS01G000030;Parent=AS01G000030

于是我尝试写了一个脚本去更改以上的错误:
`python transcript_rename.py`
脚本:
结果:
> chr1A   EVM     gene    573979  576937  .       +       .       ID=AS01G000030;Name=AS01G000030
> chr1A   EVM     mRNA    573979  576937  .       +       .       ID=AS01G000030.1;Parent=AS01G000030;Name=AS01G000030.1
> chr1A   EVM     exon    573979  574093  .       +       .       ID=AS01G000030.1.exon1;Parent=AS01G000030.1
> chr1A   EVM     CDS     573979  574093  .       +       0       ID=AS01G000030.1.cds1;Parent=AS01G000030.1
> chr1A   EVM     exon    574176  574270  .       +       .       ID=AS01G000030.1.exon2;Parent=AS01G000030.1
> chr1A   EVM     CDS     574176  574270  .       +       2       ID=AS01G000030.1.cds2;Parent=AS01G000030.1
> chr1A   EVM     exon    574381  574472  .       +       .       ID=AS01G000030.1.exon3;Parent=AS01G000030.1
> chr1A   EVM     CDS     574381  574472  .       +       0       ID=AS01G000030.1.cds3;Parent=AS01G000030.1
> chr1A   EVM     exon    574969  575062  .       +       .       ID=AS01G000030.1.exon4;Parent=AS01G000030.1
> chr1A   EVM     CDS     574969  575062  .       +       1       ID=AS01G000030.1.cds4;Parent=AS01G000030.1

