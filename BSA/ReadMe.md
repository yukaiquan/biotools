快速滑窗计算delta SNP

`python snpIndex.py chr_dens.txt chrlength.txt 1000000 output.txt`

参数依次是：

1. chr_dens.txt 密度文件 tab符分割列

chr1A   1407208 0.333333333333333
chr1A   2145558 -0.245322245322245
chr1A   2563655 0
chr1A   2563662 0.115079365079365
chr1A   3261710 0
chr1A   3278806 0
chr1A   3506660 0.230769230769231
chr1A   3507181 0

	2. chrlength.txt 染色体文件

chr1A   500000000
chr2A   460000000

 	2. 1000000 窗口大小
 	3. output.txt 输出文件