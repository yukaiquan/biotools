# JCVI多基因组共线性绘图
##### 🙈: yukaiquan
##### 📧: 1962568272@qq.com
## 1.环境配置
### 1.1 conda环境配置

```shell
conda create -n jcvi
#换中科大源
conda install -c bioconda jcvi #这个不行
conda install -y -n jcvi -c bioconda last
pip install latex
#容易报错位置
sudo apt-get install dvipng
sudo apt-get install -y texlive texlive-latex-extra texlive-latex-recommended texlive-fonts-recommended
```

## 2.准备数据

```shell
python -m jcvi.formats.gff bed --type=mRNA --key=ID ACD.sativa.gene.gff3.gz -o sfs.bed
python -m jcvi.formats.gff bed --type=mRNA --key=ID Asativa_sang.v1.1.gff3 -o sang.bed
python -m jcvi.formats.gff bed --type=mRNA --key=ID PepsiCo_OT3098_v2_gene_annotations.gff3 -o otv2.bed
python -m jcvi.formats.bed uniq sang.bed
python -m jcvi.formats.bed uniq otv2.bed
python -m jcvi.formats.bed uniq sfs.bed
seqkit grep -f <(cut -f 4 sang.uniq.bed ) sang.cds | seqkit seq -i > sang.uniq.cds
seqkit grep -f <(cut -f 4 otv2.uniq.bed ) otv2.cds | seqkit seq -i > otv2.uniq.cds
seqkit grep -f <(cut -f 4 sang.uniq.bed ) sang.pep | seqkit seq -i > sang.uniq.pep
seqkit grep -f <(cut -f 4 otv2.uniq.bed ) otv2.pep | seqkit seq -i > otv2.uniq.pep
ln -s ../otv2.uniq.bed ./otv2.bed
ln -s ../otv2.uniq.cds ./otv2.cds
ln -s ../otv2.uniq.pep ./otv2.pep
ln -s ../sang.uniq.pep ./sang.pep
ln -s ../sang.uniq.bed ./sang.bed
ln -s ../sang.uniq.cds ./sang.cds
```

## 3. 绘图

### 3.1 新建文件夹防止报错，报错清空

```shell
mkdir 01_ACD
cd 01_ACD/
ln -s ../*.bed ./
ln -s ../*.cds ./
ln -s ../*.pep ./
```

### 3.2 共线性比对绘图

```shell
python -m jcvi.compara.catalog ortholog otv2 sang --cscore=.99 --no_strip_names
python -m jcvi.compara.catalog ortholog sang sfs --cscore=.99 --no_strip_names
# build .simpe
python -m jcvi.compara.synteny screen --minspan=30 --simple otv2.sang.anchors otv2.sang.anchors.new
python -m jcvi.compara.synteny screen --minspan=30 --simple sang.sfs.anchors sang.sfs.anchors.new
python -m jcvi.graphics.karyotype seqids layout
python -m jcvi.compara.catalog ortholog otv2 sfs --cscore=.99 --no_strip_names
python -m jcvi.compara.catalog ortholog sfs sang --cscore=.99 --no_strip_names
python -m jcvi.compara.synteny screen --minspan=30 --simple otv2.sfs.anchors otv2.sfs.anchors.new
python -m jcvi.compara.synteny screen --minspan=30 --simple sfs.sang.anchors sfs.sang.anchors.new
python -m jcvi.graphics.karyotype seqids layout3
```

## 3 两基因组绘图
### 3.1 绘图代码
```bash
python -m jcvi.formats.gff bed --type=mRNA --key=ID /mnt/d/databasezip/genome/A/longiglumis/SCAU_lon/genome20220304/A_longiglumis_CN58139.gff3.gz -o scau_lon.bed
python -m jcvi.formats.gff bed --type=mRNA --key=ID /mnt/d/databasezip/genome/A/longiglumis/Sang_lon/sang_longiglumis.gff.gz -o sang_lon.bed
python -m jcvi.formats.bed uniq sang_lon.bed
python -m jcvi.formats.bed uniq scau_lon.bed

seqkit grep -f <(cut -f 4 scau_lon.uniq.bed ) A_longiglumis_CN58139_cds.fasta | seqkit seq -i > scau_lon.cds
seqkit grep -f <(cut -f 4 sang_lon.uniq.bed ) A_longiglumis_CN58138_cds.fasta | seqkit seq -i > sang_lon.cds
seqkit grep -f <(cut -f 4 scau_lon.uniq.bed ) A_longiglumis_CN58139_pep.fasta | seqkit seq -i > scau_lon.pep
seqkit grep -f <(cut -f 4 sang_lon.uniq.bed ) A_longiglumis_CN58138_pep.fasta | seqkit seq -i > sang_lon.pep


python -m jcvi.compara.catalog ortholog scau_lon sang_lon --cscore=.99 --no_strip_names
python -m jcvi.compara.synteny screen --minspan=30 --simple scau_lon.sang_lon.anchors scau_lon.sang_lon.anchors.new
python -m jcvi.graphics.karyotype seqids layout
```
### 3.2 绘图配置
绘图配置文件
seqid文件
:::danger
chr1A,chr2A,chr3A,chr4A,chr5A,chr6A,chr7A
chr1A,chr2A,chr3A,chr4A,chr5A,chr6A,chr7A
:::
layout文件
上半部分：
第一列代表在纵坐标上的位置
第二列代表横坐标的起始位置
第三列代表横坐标的终止位置
第四列代表以左边为原点的角度，向上为正向下为负
第五列标签颜色
第六列标签名
第七列位置
第八列上面生成的bed文件
下半部分：
第一列我目前还没搞明白
第二列纸袋上半部分的位置标号比如这个上半部分有两个就是分别为0 1 但是scau_lon.sang_lon.anchors.simple文件第一个指的是上半部的第一个 scau_lon 那么代号就是0

:::info
# y, xstart, xend, rotation, color, label, va,  bed
 .6,     .2,    .8,       0,      , scau_lon, top, scau_lon.bed
 .4,     .2,    .8,       0,      , sang_lon, top, sang_lon.bed
# edges
e, 0, 1, scau_lon.sang_lon.anchors.simple
:::
### 3.3 绘图
`python -m jcvi.graphics.karyotype seqids layout`

## 4 三基因组绘图
```bash
python -m jcvi.formats.gff bed --type=mRNA --key=ID /mnt/d/databasezip/genome/A/strigosa/genome20220304/A.strigosa.gene.gff3 -o strigosa.bed
python -m jcvi.formats.gff bed --type=mRNA --key=ID /mnt/d/databasezip/genome/A/atlantica/genome20220304/A_atlantica_bmc_sorted.gff3.gz -o atlantica.bed


seqkit grep -f <(cut -f 4 strigosa.uniq.bed ) A.strigosa.chr.HCC.cds | seqkit seq -i > strigosa.cds
seqkit grep -f <(cut -f 4 atlantica.uniq.bed ) Avena_atlantica.CDS.fasta | seqkit seq -i > atlantica.cds


python -m jcvi.compara.catalog ortholog atlantica strigosa --cscore=.99 --no_strip_names
python -m jcvi.compara.catalog ortholog strigosa sang_lon --cscore=.99 --no_strip_names
python -m jcvi.compara.catalog ortholog strigosa scau_lon --cscore=.99 --no_strip_names


python -m jcvi.compara.synteny screen --minspan=30 --simple strigosa.atlantica.anchors strigosa.atlantica.anchors.new
python -m jcvi.compara.synteny screen --minspan=30 --simple atlantica.scau_lon.anchors atlantica.scau_lon.anchors.new
python -m jcvi.compara.synteny screen --minspan=30 --simple atlantica.sang_lon.anchors atlantica.sang_lon.anchors.new



```
### 4.1 配置文件
**这里面才是精髓**
`seqid`
:::info
chr1A,chr2A,chr3A,chr4A,chr5A,chr6A,chr7A
chr1A,chr2A,chr3A,chr4A,chr5A,chr6A,chr7A
chr1A,chr2A,chr3A,chr4A,chr5A,chr6A,chr7A
chr1A,chr2A,chr3A,chr4A,chr5A,chr6A,chr7A
:::
`layout`
:::info
# y,	xstart,	xend,	rotation,	color,	label,	va,	bed
 .8,     .3,    .7,      0,	,	strigosa, top, strigosa.bed
 .6,     .3,    .7,       0,	,	atlantica, top, atlantica.bed
 .2,     .1,    .4,     0,	,	scau_lon, top, scau_lon.bed 
 .2,     .6,    .9,     0,	,	sang_lon, top, sang_lon.bed 
# edges
e,	0,	1,	strigosa.atlantica.anchors.simple
e,	1,	2,	atlantica.scau_lon.anchors.simple
e,	1,	3,	atlantica.sang_lon.anchors.simple
:::
在`atlantica.sang_lon.anchors.simple`在块前面加r*、g*、b*等改变颜色
绘制 `python -m jcvi.graphics.karyotype seqids layout`
### 4.2 结果
![image.png](https://cdn.nlark.com/yuque/0/2022/png/32371416/1660482561267-d5982791-c00b-4b7a-a537-5ba100a40485.png#averageHue=%23f4f3f3&clientId=u93de284a-edce-4&from=paste&height=591&id=u3e0a7a2c&name=image.png&originHeight=650&originWidth=1022&originalType=binary&ratio=1&rotation=0&showTitle=false&size=180050&status=done&style=none&taskId=u2809ddc9-e3fd-40f9-9765-3e05d9924db&title=&width=929.0908889534064)
其实还可以更骚~~~~

## 4 多个基因组共线性绘大图
上面有的数据就复用
```shell
python -m jcvi.formats.gff bed --type=mRNA --key=ID /mnt/d/databasezip/genome/ACD/SFS/07_generename/ACD_sativa_sfs.gff3.gz -o sfs.bed
python -m jcvi.formats.gff bed --type=mRNA --key=ID /mnt/d/databasezip/genome/A/atlantica/genome20220917/A_atlantica_bmc_sorted.gff3.gz -o atlan.bed
python -m jcvi.formats.gff bed --type=mRNA --key=ID /mnt/d/databasezip/genome/CD/SAU_ins/genome20220303/CD_insularis_CN108634.gff3.gz -o scauins.bed



python -m jcvi.formats.bed uniq sfs.bed
python -m jcvi.formats.bed uniq atlan.bed
python -m jcvi.formats.bed uniq scauins.bed


mv sfs.uniq.bed sfs.bed
mv atlan.uniq.bed atlan.bed
mv scauins.uniq.bed scauins.bed


seqkit grep -f <(cut -f 4 sfs.bed ) /mnt/d/databasezip/genome/ACD/SFS/07_generename/ACD_sativa_sfs_cds.fasta.gz | seqkit seq -i > sfs.cds
seqkit grep -f <(cut -f 4 atlan.bed ) /mnt/d/databasezip/genome/A/atlantica/genome20220917/A_atlantica_bmc_cds.fasta.gz | seqkit seq -i > atlan.cds
seqkit grep -f <(cut -f 4 scauins.bed ) /mnt/d/databasezip/genome/CD/SAU_ins/genome20220303/CD_insularis_CN108634_cds.fasta.gz | seqkit seq -i > scauins.cds

seqkit grep -f <(cut -f 4 atlan.bed ) /mnt/d/databasezip/genome/A/atlantica/genome20220917/A_atlantica_bmc_pep.fasta.gz | seqkit seq -i > atlan.pep
seqkit grep -f <(cut -f 4 sfs.bed ) /mnt/d/databasezip/genome/ACD/SFS/07_generename/ACD_sativa_sfs_pep.fasta.gz | seqkit seq -i > sfs.pep
seqkit grep -f <(cut -f 4 scauins.bed ) /mnt/d/databasezip/genome/CD/SAU_ins/genome20220303/CD_insularis_CN108634_pep.fasta.gz | seqkit seq -i > scauins.pep


##上面这些代码太繁琐了作为一个懒狗 得封装一下##
python /mnt/e/software/Ex-seq/biotools/jcvi/jcvi_files_product.py /mnt/d/databasezip/genome/CD/Sang_ins/CD_insularis_BYU209.gff3.gz /mnt/d/databasezip/genome/CD/Sang_ins/CD_insularis_BYU209_cds.fasta.gz /mnt/d/databasezip/genome/CD/Sang_ins/CD_insularis_BYU209_pep.fasta.gz sangins
python /mnt/e/software/Ex-seq/biotools/jcvi/jcvi_files_product.py /mnt/d/databasezip/genome/A/strigosa/genome20220917/A_strigosa_nc_sorted.gff3.gz /mnt/d/databasezip/genome/A/strigosa/genome20220917/A_strigosa_nc_cds.fasta.gz /mnt/d/databasezip/genome/A/strigosa/genome20220917/A_strigosa_nc_pep.fasta.gz ncstri
python /mnt/e/software/Ex-seq/biotools/jcvi/jcvi_files_product.py /mnt/d/databasezip/genome/A/eriantha/genome_rename/C_eriantha_bmc.gff3.gz /mnt/d/databasezip/genome/A/eriantha/genome_rename/C_eriantha_bmc_cds.fasta.gz /mnt/d/databasezip/genome/A/eriantha/genome_rename/C_eriantha_bmc_pep.fasta.gz bmceri
python /mnt/e/software/Ex-seq/biotools/jcvi/jcvi_files_product.py /mnt/d/databasezip/genome/A/longiglumis/SCAU_lon/genome20220304/A_longiglumis_CN58139.gff3.gz /mnt/d/databasezip/genome/A/longiglumis/SCAU_lon/genome20220304/A_longiglumis_CN58139_cds.fasta.gz /mnt/d/databasezip/genome/A/longiglumis/SCAU_lon/genome20220304/A_longiglumis_CN58139_pep.fasta.gz scaulon
python /mnt/e/software/Ex-seq/biotools/jcvi/jcvi_files_product.py /mnt/d/databasezip/genome/A/longiglumis/Sang_lon/A_longiglumis_CN58138.gff3.gz /mnt/d/databasezip/genome/A/longiglumis/Sang_lon/A_longiglumis_CN58138_cds.fasta.gz /mnt/d/databasezip/genome/A/longiglumis/Sang_lon/A_longiglumis_CN58138_pep.fasta.gz sanglon


# 基于CDS共线性
python -m jcvi.compara.catalog ortholog otv2 sang --cscore=.99 --no_strip_names
python -m jcvi.compara.catalog ortholog sang sfs --cscore=.99 --no_strip_names

python -m jcvi.compara.catalog ortholog sfs sangins --cscore=.99 --no_strip_names
python -m jcvi.compara.catalog ortholog sfs scauins --cscore=.99 --no_strip_names

python -m jcvi.compara.catalog ortholog sfs scaulon --cscore=.99 --no_strip_names
# 需要构思一下最终的图
python -m jcvi.compara.catalog ortholog sfs sanglon --cscore=.99 --no_strip_names
python -m jcvi.compara.catalog ortholog sanglon ncstri --cscore=.99 --no_strip_names
python -m jcvi.compara.catalog ortholog ncstri atlan --cscore=.99 --no_strip_names
python -m jcvi.compara.catalog ortholog scaulon sanglon --cscore=.99 --no_strip_names
python -m jcvi.compara.catalog ortholog sangins bmceri --cscore=.99 --no_strip_names
python -m jcvi.compara.catalog ortholog scauins sangins --cscore=.99 --no_strip_names


python -m jcvi.compara.synteny screen --minspan=30 --simple otv2.sang.anchors otv2.sang.anchors.new
python -m jcvi.compara.synteny screen --minspan=30 --simple sang.sfs.anchors sang.sfs.anchors.new
python -m jcvi.compara.synteny screen --minspan=30 --simple sfs.scauins.anchors sfs.scauins.anchors.new
python -m jcvi.compara.synteny screen --minspan=30 --simple scauins.sangins.anchors scauins.sangins.anchors.new
python -m jcvi.compara.synteny screen --minspan=30 --simple sangins.bmceri.anchors sangins.bmceri.anchors.new
python -m jcvi.compara.synteny screen --minspan=30 --simple sfs.sanglon.anchors sfs.sanglon.anchors.new
python -m jcvi.compara.synteny screen --minspan=30 --simple sfs.scaulon.anchors sfs.scaulon.anchors.new
python -m jcvi.compara.synteny screen --minspan=30 --simple scaulon.sanglon.anchors scaulon.sanglon.anchors.new
python -m jcvi.compara.synteny screen --minspan=30 --simple sanglon.ncstri.anchors sanglon.ncstri.anchors.new
python -m jcvi.compara.synteny screen --minspan=30 --simple ncstri.atlan.anchors ncstri.atlan.anchors.new

```
`python -m jcvi.graphics.karyotype seqids layout`
seqids
> chr1A,chr2A,chr3A,chr4A,chr5A,chr6A,chr7A,chr1C,chr2C,chr3C,chr4C,chr5C,chr6C,chr7C,chr1D,chr2D,chr3D,chr4D,chr5D,chr6D,chr7D
> chr1A,chr2A,chr3A,chr4A,chr5A,chr6A,chr7A,chr1C,chr2C,chr3C,chr4C,chr5C,chr6C,chr7C,chr1D,chr2D,chr3D,chr4D,chr5D,chr6D,chr7D
> chr1A,chr2A,chr3A,chr4A,chr5A,chr6A,chr7A,chr1C,chr2C,chr3C,chr4C,chr5C,chr6C,chr7C,chr1D,chr2D,chr3D,chr4D,chr5D,chr6D,chr7D
> chr1C,chr2C,chr3C,chr4C,chr5C,chr6C,chr7C,chr1D,chr2D,chr3D,chr4D,chr5D,chr6D,chr7D
> chr1C,chr2C,chr3C,chr4C,chr5C,chr6C,chr7C,chr1D,chr2D,chr3D,chr4D,chr5D,chr6D,chr7D
> chr1C,chr2C,chr3C,chr4C,chr5C,chr6C,chr7C
> chr1A,chr2A,chr3A,chr4A,chr5A,chr6A,chr7A
> chr1A,chr2A,chr3A,chr4A,chr5A,chr6A,chr7A
> chr1A,chr2A,chr3A,chr4A,chr5A,chr6A,chr7A
> chr1A,chr2A,chr3A,chr4A,chr5A,chr6A,chr7A

layout：
> # y,    xstart, xend,   rotation,       color,  label,  va,     bed
>  .9,   .2, .8,      0,	,       otv2, top, otv2.bed
>  .8,   .2, .8,       0,        ,       sang, top, sang.bed
>  .7,   .2, .8,     0,	,       sfs, top,   sfs.bed
>  .5,   .5, .9, 0,	,   scauins,    top,    scauins.bed
>  .3,   .5, .9, 0,	,   sangins,    top,    sangins.bed
>  .1,   .5, .7, 0,	,	bmceri,	top,	bmceri.bed
>  .4,	.1,	.3,	0,	,	scaulon,	top,	scaulon.bed
>  .3,	.1,	.3,	0,	,	sanglon,	top,	sanglon.bed	
>  .2,	.1,	.3,	0,	,	ncstri,	top,	ncstri.bed
>  .1,	.1,	.3,	0,	,	atlan,	top,	atlan.bed
> # edges
> e,      0,      1,      otv2.sang.anchors.simple
> e,      1,      2,      sang.sfs.anchors.simple
> e,      2,      3,      sfs.scauins.anchors.simple
> e,      3,      4,      scauins.sangins.anchors.simple
> e,      4,      5,      sangins.bmceri.anchors.simple
> e,      2,      6,      sfs.scaulon.anchors.simple
> e,      6,      7,      scaulon.sanglon.anchors.simple
> e,      7,      8,      sanglon.ncstri.anchors.simple
> e,      8,      9,      ncstri.atlan.anchors.simple

![oatgenomeall.png](https://cdn.nlark.com/yuque/0/2023/png/32371416/1676772851734-9092e92b-1598-440f-a11f-bc4dc97d0150.png#averageHue=%23efefef&clientId=ud632c3de-0eb6-4&from=drop&id=ub85b4623&name=oatgenomeall.png&originHeight=1400&originWidth=1600&originalType=binary&ratio=1&rotation=0&showTitle=false&size=99785&status=done&style=none&taskId=u0396e416-53fd-4ce7-b223-ce5c2a09906&title=)
