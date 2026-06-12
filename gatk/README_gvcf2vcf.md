# GVCF to VCF Pipeline

分区域合并GVCF并生成VCF文件的流水线工具。

## 功能概述

该工具分为三个主要步骤：

### Step 0: 生成基因组区间文件

从基因组的 `.fai` 文件生成区间文件，用于后续的分区域处理。

**区间文件格式（带重叠）：**

```
chr1H  1        5000000   1        4990000
chr1H  4980001  10000000  4990001  9990000
chr1H  9980001  15000000  9990001  14990000
...
```

- 第1列：染色体名称
- 第2-3列：**提取区间**（用于 bcftools view -r，带重叠）
- 第4-5列：**VCF过滤区间**（非重叠，用于最终合并）

**重叠区间原理：**

- **提取区间重叠**：4980001-10000000 与前一区间重叠 20000 bp，确保边界变异被完整捕获
- **过滤区间不重叠**：4990001-9990000 紧接前一个区间结束位置，避免重复记录

**示例：**
对于 5Mb 窗口和 20kb 重叠：

```
窗口1: 提取 chr1H:1-5000000,        过滤 chr1H:1-4990000
窗口2: 提取 chr1H:4980001-10000000, 过滤 chr1H:4990001-9990000
窗口3: 提取 chr1H:9980001-15000000, 过滤 chr1H:9990001-14990000
```

### Step 1: 分区域提取并合并GVCF

对每个区间：

1. 使用 bcftools 从每个样本的 GVCF 中提取该区域的变异
2. 使用 GATK CombineGVCFs 合并所有样本的区域 GVCF
3. 清理中间文件

### Step 2: 基因分型并生成VCF

1. 合并所有区间的 GVCF 文件
2. 使用 GATK GenotypeGVCFs 进行基因分型
3. 分离 SNP 和 INDEL
4. 按染色体合并区间 VCF 文件
5. 清理中间文件

## 使用方法

### 方式一：逐步执行

#### Step 0: 生成区间文件

```bash
python gvcf2vcf_tools.py step0 \
    -f genome.fa.fai \
    -o intervals.txt \
    -w 5000000 \
    --overlap 20000
```

**参数说明：**

- `-f, --fai`: 基因组的 fai 索引文件
- `-o, --output`: 输出的区间文件路径
- `-w, --window`: 窗口大小（默认：5000000）
- `--overlap`: 区间重叠大小（默认：20000）
- `-o, --output`: 输出的区间文件路径
- `-w, --window`: 窗口大小（默认：5000000）

#### Step 1: 提取并合并GVCF

```bash
python gvcf2vcf_tools.py step1 \
    -i gvcf_list.txt \
    -l intervals.txt \
    -o interval_gvcfs/ \
    -r reference.fa \
    -t 8 \
    --java-options "-Xmx40g"
```

**参数说明：**

- `-i, --input`: GVCF 文件列表
- `-l, --intervals`: Step 0 生成的区间文件
- `-o, --output`: 输出目录
- `-r, --reference`: 参考基因组文件
- `-t, --thread`: 线程数（默认：8）
- `--java-options`: GATK Java 参数（默认：-Xmx40g）
- `--dry-run`: 仅打印命令不执行
- `--keep-intermediate`: 保留中间文件

**GVCF 列表文件格式：**

```
HOR_16340    ../../../s1.gvcf/HOR_16340/gvcf.gz
HOR_9582     ../../../s1.gvcf/HOR_9582/gvcf.gz
...
```

或者每行只包含路径：

```
../../../s1.gvcf/HOR_16340/gvcf.gz
../../../s1.gvcf/HOR_9582/gvcf.gz
...
```

#### Step 2: 基因分型

```bash
python gvcf2vcf_tools.py step2 \
    -d interval_gvcfs/ \
    -l intervals.txt \
    -o final_vcfs/ \
    -r reference.fa \
    -t 8 \
    --java-options "-Xmx40g" \
    -p cohort
```

**参数说明：**

- `-d, --dir`: Step 1 生成的区间 GVCF 目录
- `-l, --intervals`: 区间文件
- `-o, --output`: 最终 VCF 输出目录
- `-r, --reference`: 参考基因组文件
- `-t, --thread`: 线程数
- `--java-options`: GATK Java 参数
- `-p, --prefix`: 输出文件前缀（默认：cohort）
- `--dry-run`: 仅打印命令不执行
- `--keep-intermediate`: 保留中间文件

### 方式二：一键执行所有步骤

```bash
python gvcf2vcf_tools.py all \
    -i gvcf_list.txt \
    -f genome.fa.fai \
    -r reference.fa \
    -o output/ \
    -t 8 \
    -w 5000000 \
    --overlap 20000 \
    --java-options "-Xmx40g" \
    -p cohort
```

**输出目录结构：**

```
output/
├── intervals.txt           # 区间文件
├── interval_gvcfs/         # Step 1 中间结果
│   ├── 1.gvcf             # 区间1合并后的GVCF
│   ├── 2.gvcf
│   └── ...
└── final_vcfs/             # Step 2 最终结果
    ├── cohort.snp.vcf.gz   # 全基因组SNP
    ├── cohort.indel.vcf.gz # 全基因组INDEL
    ├── chr1H.snp.vcf.gz    # 按染色体分组的SNP
    ├── chr1H.indel.vcf.gz
    └── ...
```

## 示例工作流程

### 示例 1: 处理大麦重测序数据

假设你有以下文件：

- 参考基因组：`/path/to/Morex.fa`
- GVCF 列表：`gvcf_list.txt`
- 基因组索引：`/path/to/Morex.fa.fai`

```bash
# 完整流程
python gvcf2vcf_tools.py all \
    -i gvcf_list.txt \
    -f /path/to/Morex.fa.fai \
    -r /path/to/Morex.fa \
    -o barley_wgs_output/ \
    -t 16 \
    -w 5000000 \
    --java-options "-Xmx80g"
```

### 示例 2: 分步调试

```bash
# Step 0: 生成区间（5Mb窗口，20kb重叠）
python gvcf2vcf_tools.py step0 \
    -f Morex.fa.fai \
    -o intervals.txt \
    -w 5000000 \
    --overlap 20000

# Step 1: 提取区域（先测试，不实际运行）
python gvcf2vcf_tools.py step1 \
    -i gvcf_list.txt \
    -l intervals.txt \
    -o interval_gvcfs/ \
    -r Morex.fa \
    --dry-run

# Step 1: 实际运行
python gvcf2vcf_tools.py step1 \
    -i gvcf_list.txt \
    -l intervals.txt \
    -o interval_gvcfs/ \
    -r Morex.fa \
    -t 16

# Step 2: 基因分型
python gvcf2vcf_tools.py step2 \
    -d interval_gvcfs/ \
    -l intervals.txt \
    -o final_vcfs/ \
    -r Morex.fa \
    -p barley_panel
```

## 核心算法说明

### 区间重叠处理

为了确保边界处的变异被完整捕获且不重复，采用以下策略：

**双重坐标系统：**

- **提取区间（带重叠）**：使用 `bcftools view -r chr:start-end` 提取
- **过滤区间（不重叠）**：在合并最终VCF时，使用 `awk '$2>=vcf_start && $2<=vcf_end'` 过滤

**示例（5Mb窗口，20kb重叠）：**

```
区间1:
  提取: chr1H:1-5000000       (完整窗口)
  过滤: chr1H:1-4990000       (去掉最后20kb)

区间2:
  提取: chr1H:4980001-10000000 (与前区间重叠20kb)
  过滤: chr1H:4990001-9990000  (紧接区间1的过滤结束位置)

区间3:
  提取: chr1H:9980001-15000000 (与前区间重叠20kb)
  过滤: chr1H:9990001-14990000 (紧接区间2的过滤结束位置)
```

**工作流程：**

1. **Step 1 提取阶段**：使用重叠的提取区间，确保边界变异（如 4990001-5000000）在相邻区间都被捕获
2. **Step 2 合并阶段**：使用非重叠的过滤区间，每个位置只保留一次，避免重复记录

**优势：**

- ✅ 边界变异不会丢失
- ✅ 最终VCF没有重复记录
- ✅ 重叠区域作为缓冲区，确保变异检测的连续性

### 并行处理建议

虽然脚本本身是串行的，但可以利用以下方式并行化：

1. **不同染色体并行**：在不同节点上处理不同染色体
2. **区间并行**：修改脚本，使不同区间可以并行处理

## 依赖工具

- **GATK** (Genome Analysis Toolkit)
  - CombineGVCFs
  - GenotypeGVCFs
  - SelectVariants
- **bcftools**: 用于提取GVCF区域

- **bgzip/htslib**: 用于压缩VCF文件

## 注意事项

1. **内存需求**：GATK CombineGVCFs 和 GenotypeGVCFs 需要大量内存，建议根据样本数量调整 `--java-options`

2. **磁盘空间**：中间文件可能占用大量空间，建议定期清理（默认自动清理）

3. **断点续传**：脚本会检查输出文件是否存在且非空，如果已存在则跳过该步骤

4. **样本数量**：大量样本（>100）建议增加内存和减小窗口大小

## 常见问题

### Q1: 如何处理不同大小的染色体？

A: 脚本会自动处理，最后一个窗口会根据染色体长度调整。

### Q2: 可以跳过某些步骤吗？

A: 可以，使用分步模式。例如，如果已有区间文件，可以直接从 Step 1 开始。

### Q3: 如何验证结果的正确性？

A: 检查最终 VCF 文件的变异位点数量和样本数是否符合预期。

## 更新日志

### v2.0 (2026-06-05)

- 完全重构为三步流水线
- 新增区间文件生成功能
- 优化内存使用和中间文件管理
- 支持断点续传
- 添加详细的日志输出

## 作者

biotools 开发团队

## 许可证

MIT License
