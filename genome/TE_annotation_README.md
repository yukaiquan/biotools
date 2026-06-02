# TE注释流程 (TE Annotation Pipeline)

## 概述

本流程用于对基因组进行转座子（Transposable Element, TE）注释，采用三步骤策略：

1. **MITE-Hunter**: 搜索MITE（Miniature Inverted-repeat Transposable Elements）转座子
2. **LTR注释**: LTR_Finder + LTR_Harvest + LTR_Retriever 构建LTR库
3. **RepeatModeler/RepeatMasker**: de novo搜索 + 全基因组重复序列注释

## 流程图

```
TR softmasked genome
       ↓
[Step 1] MITE-Hunter → MITE库 (MITE.fa)
       ↓
[Step 2] LTR_Finder (并行) ┐
       ↓                  ├→ LTR_Retriever → LTR库 (LTR.fa)
       LTR_Harvest ────────┘
       ↓
[Step 3] MITE.fa + LTR.fa → TE.lib → Hardmask → RepeatModeler → RepMod.lib
       ↓                                                      ↓
       TEclass分类Unknown ←───────────────────────────────────┘
       ↓
       TE.lib + RepMod.lib + Repbase → 最终TE库 (final.lib)
       ↓
       RepeatMasker全基因组注释 → 最终结果
```

## 软件依赖

### 必需软件

| 软件              | 版本          | 用途                  | 参考文献                   |
| ----------------- | ------------- | --------------------- | -------------------------- |
| **MITE-Hunter**   | -             | 搜索MITE转座子        | Han Y et al; 2010          |
| **LTR_Finder**    | -             | LTR反转录转座子识别   | Xu Z et al; 2007           |
| **LTR_Harvest**   | (GenomeTools) | LTR反转录转座子识别   | D. Ellinghaus et al; 2008  |
| **LTR_Retriever** | -             | 整合LTR结果，去假阳性 | Ou S et al; 2018           |
| **RepeatModeler** | -             | de novo重复序列搜索   | Bedell JA et al; 2000      |
| **TEclass**       | -             | TE分类工具            | György Abrusán et al; 2009 |
| **RepeatMasker**  | -             | 重复序列注释          | Smit AFA et al; 2013       |

### 辅助工具

- `samtools`: 基因组索引和序列提取
- `cd-hit-est`: 序列去重
- `seqtk`: 序列处理
- `gt` (GenomeTools): LTR_Harvest依赖

### 可选数据库

- **Repbase**: 已知重复序列数据库 (J, J. et al; 2005)

## 安装配置

### 1. 修改软件路径

编辑脚本顶部的软件路径配置：

```python
# MITE-Hunter
MITEHUNTER = "/path/to/MITE-Hunter/MITE_Hunter_manager.pl"

# LTR相关工具
LTR_FINDER = "/path/to/ltr_finder"
LTR_HARVEST = "gt"  # GenomeTools
LTR_RETRIEVER = "/path/to/LTR_retriever/LTR_retriever.pl"

# RepeatModeler
REPEATMODELER = "/path/to/RepeatModeler/RepeatModeler"
BUILD_DATABASE = "/path/to/RepeatModeler/BuildDatabase"

# TEclass
TECLASS = "/path/to/TEclass/TEclassTest.pl"

# RepeatMasker
REPEATMASKER = "/path/to/RepeatMasker"
REPEATMASKER_LIB = "/path/to/RepBase/RepBase.lib"  # 可选
```

### 2. 安装依赖

```bash
# Conda安装（推荐）
conda install -c bioconda ltr_finder cd-hit repeatmasker RepeatModeler

# 手动安装MITE-Hunter
wget http://target.iplantcollaborative.org/mitehunter.html
tar -xzf MITE-Hunter.tar.gz
cd MITE-Hunter
# 配置Perl依赖

# 安装LTR_Retriever
git clone https://github.com/oushujun/LTR_retriever.git
cd LTR_retriever
# 配置路径

# 安装TEclass
git clone https://github.com/labfei/TEclass.git
cd TEclass
# 配置Python/R依赖
```

## 使用方法

### 基本用法

```bash
# 运行完整流程（所有3个步骤）
python TE_annotation.py -g genome.fasta -s species_name -t 16

# 参数说明：
#   -g, --genome: 输入基因组文件（TR softmasked）
#   -s, --species: 物种名称（用于输出文件命名）
#   -t, --threads: 线程数（默认: 8）
```

### 分步运行

#### Step 1: MITE-Hunter

```bash
# 只运行Step 1
python TE_annotation.py -g genome.fasta -s species_name --step1

# 输出文件：
#   01_MITE/species_name.MITE.fa
```

**采样比例计算**：

- 基因组 ≤ 700 Mb: 使用全部序列 (ratio = 1.0)
- 基因组 > 700 Mb: ratio = 1 / (genome_size / 700)
- 示例: 3 Gb基因组 → ratio ≈ 0.23

#### Step 2: LTR注释

```bash
# 只运行Step 2
python TE_annotation.py -g genome.fasta -s species_name --step2

# 输出文件：
#   02_LTR_Finder/species_name.ltr.finder_*.scn
#   03_LTR_Harvest/species_name.harvest.scn
#   04_LTR_Retriever/species_name.LTR.lib.fa
```

**并行策略**：

- LTR_Finder按染色体并行运行
- 每个染色体使用指定线程数

#### Step 3: RepeatModeler + RepeatMasker

```bash
# 只运行Step 3（需要提供前两步的库文件）
python TE_annotation.py -g genome.fasta -s species_name --step3 \
    --mite-lib 01_MITE/species_name.MITE.fa \
    --ltr-lib 04_LTR_Retriever/species_name.LTR.lib.fa

# 输出文件：
#   08_TE_Library/species_name.TE.lib
#   05_RepeatModeler/species_name.RepMod.lib
#   08_TE_Library/species_name.final.lib
#   07_RepeatMasker/ (RepeatMasker结果)
```

**流程说明**：

1. 整合MITE和LTR库 → TE.lib
2. 使用TE.lib对基因组hardmask
3. RepeatModeler在hardmask基因组上de novo搜索
4. TEclass对Unknown序列分类
5. 整合TE.lib + RepMod.lib + Repbase → 最终库
6. RepeatMasker全基因组注释

## 输出文件说明

### 目录结构

```
.
├── 01_MITE/                          # Step 1输出
│   └── species_name.MITE.fa          # MITE库文件
├── 02_LTR_Finder/                    # Step 2.1输出
│   ├── species_name.ssrMasked_chr*.fasta
│   └── species_name.ltr.finder_*.scn
├── 03_LTR_Harvest/                   # Step 2.2输出
│   └── species_name.harvest.scn
├── 04_LTR_Retriever/                 # Step 2.3输出
│   └── species_name.LTR.lib.fa       # LTR库文件
├── 05_RepeatModeler/                 # Step 3.3输出
│   └── species_name.RepMod.lib       # RepeatModeler库
├── 06_TEclass/                       # Step 3.4输出
│   └── species_name.classified.lib
├── 07_RepeatMasker/                  # Step 3.6输出
│   ├── genome.fasta.masked           # softmask基因组
│   ├── genome.fasta.out              # RepeatMasker表格结果
│   └── genome.fasta.out.gff          # GFF格式结果
├── 08_TE_Library/                    # TE库文件
│   ├── species_name.TE.lib           # MITE + LTR整合库
│   └── species_name.final.lib        # 最终库（含Repbase）
└── species_name.TE_annotation.log    # 日志文件
```

### 主要结果文件

| 文件                      | 说明                         |
| ------------------------- | ---------------------------- |
| `species_name.MITE.fa`    | MITE转座子序列库             |
| `species_name.LTR.lib.fa` | LTR反转录转座子序列库        |
| `species_name.TE.lib`     | MITE + LTR整合库             |
| `species_name.RepMod.lib` | RepeatModeler de novo库      |
| `species_name.final.lib`  | **最终TE库**（用于下游分析） |
| `genome.fasta.out`        | RepeatMasker表格结果         |
| `genome.fasta.out.gff`    | GFF格式TE注释                |

## 参数调优

### MITE-Hunter参数

```python
# 调整采样比例（脚本自动计算）
-P 0.22  # 对于3Gb基因组

# 调整候选MITE数量
-n 20    # 每个家族最多20个候选
```

### LTR_Finder参数

```python
# LTR长度范围
-L 7000  # 最大LTR长度
-l 100   # 最小LTR长度

# TSD长度
-D 15000 # 最大距离
-d 1000  # 最小距离
```

### LTR_Harvest参数

```python
# 相似度阈值
-similar 90  # LTR对相似度≥90%

# TSD和Motif
-mintsd 4    # 最小TSD长度
-maxtsd 6    # 最大TSD长度
-motif TGCA  # LTR motif
```

### RepeatMasker参数

```python
# 输出格式
-gff     # 输出GFF格式
-xsmall  # softmask（保留小写）

# 注释模式
-species species_name  # 物种特异性注释
```

## 常见问题

### 1. MITE-Hunter运行时间过长

**原因**: 基因组太大，采样比例过低

**解决**:

- 增加线程数 `-t 16`
- 手动调整采样比例（编辑脚本中的 `calculate_sample_ratio()`）

### 2. LTR_Finder内存不足

**原因**: 单个染色体过大

**解决**:

- 减少每个LTR_Finder进程的线程数
- 分批处理染色体

### 3. RepeatModeler找不到输出文件

**原因**: RepeatModeler版本差异，输出文件名不一致

**解决**:

- 检查RepeatModeler输出目录
- 手动指定输出文件（编辑脚本中的 `possible_outputs`）

### 4. TEclass分类失败

**原因**: TEclass依赖的Python/R环境问题

**解决**:

- 检查TEclass依赖
- 跳过TEclass步骤，使用未分类的RepMod.lib

### 5. RepeatMasker库文件格式错误

**原因**: 库文件包含非标准字符或格式

**解决**:

- 使用 `cd-hit-est` 去重和清理
- 手动检查库文件格式

## 引用文献

如果您使用本流程，请引用以下文献：

### MITE-Hunter

Han Y, Qin P, Weng S, et al. **MITE-Hunter: a computational approach to identify miniature inverted repeat transposable elements in plant genomes**[J]. _Bioinformatics_, 2010, 26(22): 2918-2920.

### LTR_Finder

Xu Z, Wang H. **LTR_FINDER: an efficient tool for the prediction of full-length LTR retrotransposons**[J]. _Nucleic Acids Research_, 2007, 35(suppl_2): W265-W268.

### LTR_Harvest

Ellinghaus D, Kurtz S, Willhoeft U. **LTRharvest, an efficient and generic software for de novo detection of LTR retrotransposons**[J]. _BMC Bioinformatics_, 2008, 9(1): 18.

### LTR_Retriever

Ou S, Chen J, Jiang N. **Assessing genome assembly quality using the LTR Assembly Index (LAI)**[J]. _Nucleic Acids Research_, 2018, 46(21): e126-e126.

### RepeatModeler

Bedell JA, Korf I, Gish W. **MaskerAid: a performance enhancement to RepeatMasker**[J]. _Bioinformatics_, 2000, 16(11): 1040-1041.

### TEclass

Abrusán G, Grundmann N, DeMester L, et al. **TEclass—a tool for automated classification of unknown eukaryotic transposable elements**[J]. _Bioinformatics_, 2009, 25(10): 1329-1330.

### Repbase

Jurka J, Kapitonov VV, Pavlicek A, et al. **Repbase Update, a database of eukaryotic repetitive elements**[J]. _Cytogenetic and Genome Research_, 2005, 110(1-4): 462-467.

### RepeatMasker

Smit AFA, Hubley R, Green P. **RepeatMasker Open-4.0**[J]. 2013-2015. <http://www.repeatmasker.org>.

## 版本历史

### v2.0 (2026-05-13)

- 重构为三步骤流程
- 添加MITE-Hunter支持
- 添加LTR_Retriever整合
- 添加RepeatModeler和TEclass支持
- 添加并行处理支持
- 改进日志和错误处理

### v1.0 (2023-07-18)

- 初始版本
- 基本TE注释功能

## 许可证

MIT License

## 联系方式

如有问题，请联系：[your-email@example.com]
