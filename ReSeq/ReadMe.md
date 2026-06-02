# NGS Minimap2 Pipeline

## Overview

NGS sequencing data alignment pipeline supporting both single-end and paired-end reads. The pipeline performs quality control, alignment, sorting, deduplication, and format conversion.

## Pipeline Steps

1. **Fastp QC** (Optional) - Quality control and adapter trimming
2. **Minimap2 Alignment** - Align reads to reference genome
3. **Samtools Sort** - Sort BAM file by coordinates
4. **Sambamba Markdup** (Optional) - Mark/remove duplicates
5. **BAM to CRAM** (Optional) - Convert to space-efficient CRAM format

## Features

- ✅ Supports single-end and paired-end sequencing data
- ✅ Optional Fastp quality control step
- ✅ Unique mapping filter option
- ✅ Optional duplicate marking
- ✅ Supports BAM or CRAM output format
- ✅ Automatic intermediate file cleanup
- ✅ Dry-run mode for testing

## Requirements

- Python 3.6+
- Fastp
- Minimap2
- Samtools
- Sambamba (for markdup)

## Installation

```bash
# Install via conda
conda install -c bioconda fastp minimap2 samtools sambamba
```

## Usage

### Basic Usage (Paired-end)

```bash
python ngs_minimap2.py \
    -i sample_R1.fq.gz \
    -j sample_R2.fq.gz \
    -o output_dir \
    -r reference.fa \
    -t 16 \
    -p sample \
    --rg-id sample1 \
    --rg-lb lib1 \
    --rg-sm sample1
```

### Single-end Mode

```bash
python ngs_minimap2.py \
    -i sample.fq.gz \
    -o output_dir \
    -r reference.fa \
    -t 16 \
    -p sample \
    --rg-id sample1 \
    --rg-lb lib1 \
    --rg-sm sample1
```

### Skip Fastp QC

```bash
python ngs_minimap2.py \
    -i sample_R1.fq.gz \
    -j sample_R2.fq.gz \
    -o output_dir \
    -r reference.fa \
    --skip-fastp \
    --rg-id sample1 \
    --rg-lb lib1 \
    --rg-sm sample1
```

### Unique Mapping Only

```bash
python ngs_minimap2.py \
    -i sample_R1.fq.gz \
    -j sample_R2.fq.gz \
    -o output_dir \
    -r reference.fa \
    --uniq-mapping \
    --rg-id sample1 \
    --rg-lb lib1 \
    --rg-sm sample1
```

## Parameters

### Required Parameters

- `-i, --input1`: R1 fastq file (or single-end input)
- `-o, --output_dir`: Output directory
- `-r, --reference`: Reference genome (FASTA format, indexed)
- `--rg-id`: Read group ID
- `--rg-lb`: Read group library
- `--rg-sm`: Read group sample name

### Optional Parameters

- `-j, --input2`: R2 fastq file (paired-end mode)
- `-t, --thread`: Number of threads (default: 8)
- `-p, --prefix`: Output prefix (default: sample)
- `--rg-pl`: Read group platform (default: ILLUMINA)
- `-f, --output-format`: Output format, cram or bam (default: cram)

### Flag Options

- `--dry-run`: Print commands without execution
- `--keep-intermediate`: Keep intermediate files
- `--no-markdup`: Skip duplicate marking step
- `--uniq-mapping`: Filter non-unique mappings (samtools view -q 20 -F 3332)
- `--skip-fastp`: Skip Fastp QC and use raw fastq for alignment

## Output Files

### Fastp Step (if not skipped)

- `{prefix}_fastp.html`: QC report in HTML
- `{prefix}_fastp.json`: QC report in JSON
- `{prefix}_clean_1.fq.gz`: Cleaned R1 reads
- `{prefix}_clean_2.fq.gz`: Cleaned R2 reads (paired-end only)

### Alignment Step

- `{prefix}.sorted.bam`: Sorted BAM file
- `{prefix}.sorted.bam.csi`: BAM index
- `{prefix}_minimap2.log`: Minimap2 log file

### Markdup Step (if enabled)

- `{prefix}.rmdup.bam`: Deduplicated BAM file
- `{prefix}.rmdup.bam.bai`: BAM index

### Final Output

- CRAM mode: `{prefix}_sort.cram` and `{prefix}_sort.cram.crai`
- BAM mode: Final BAM file with index

## Examples

### Standard Pipeline with All Steps

```bash
python ngs_minimap2.py \
    -i R1.fq.gz \
    -j R2.fq.gz \
    -o results \
    -r genome.fa \
    -t 32 \
    -p my_sample \
    --rg-id my_sample \
    --rg-lb my_library \
    --rg-sm my_sample \
    --output-format cram
```

### Quick Alignment (No QC, No Dedup, BAM Output)

```bash
python ngs_minimap2.py \
    -i R1.fq.gz \
    -j R2.fq.gz \
    -o results \
    -r genome.fa \
    --skip-fastp \
    --no-markdup \
    --output-format bam \
    --rg-id my_sample \
    --rg-lb my_library \
    --rg-sm my_sample
```

## Notes

1. Reference genome must be indexed (`samtools faidx reference.fa`)
2. Sambamba is required for markdup functionality
3. Intermediate files are automatically cleaned unless `--keep-intermediate` is specified
4. Use `--dry-run` to preview commands before execution

---

# NGS Minimap2 流程

## 概述

NGS 测序数据比对流程，支持单端和双端测序数据。流程包含质量控制、比对、排序、去重和格式转换等步骤。

## 流程步骤

1. **Fastp质控**（可选）- 质量控制和接头去除
2. **Minimap2比对** - 将reads比对到参考基因组
3. **Samtools排序** - 按坐标排序BAM文件
4. **Sambamba去重**（可选）- 标记/去除重复序列
5. **BAM转CRAM**（可选）- 转换为节省空间的CRAM格式

## 特性

- ✅ 支持单端和双端测序数据
- ✅ 可选的Fastp质量控制步骤
- ✅ 唯一比对过滤选项
- ✅ 可选的重复序列标记
- ✅ 支持BAM或CRAM输出格式
- ✅ 自动清理中间文件
- ✅ 测试模式的dry-run功能

## 环境要求

- Python 3.6+
- Fastp
- Minimap2
- Samtools
- Sambamba（用于去重）

## 安装

```bash
# 使用conda安装
conda install -c bioconda fastp minimap2 samtools sambamba
```

## 使用方法

### 基本用法（双端）

```bash
python ngs_minimap2.py \
    -i sample_R1.fq.gz \
    -j sample_R2.fq.gz \
    -o output_dir \
    -r reference.fa \
    -t 16 \
    -p sample \
    --rg-id sample1 \
    --rg-lb lib1 \
    --rg-sm sample1
```

### 单端模式

```bash
python ngs_minimap2.py \
    -i sample.fq.gz \
    -o output_dir \
    -r reference.fa \
    -t 16 \
    -p sample \
    --rg-id sample1 \
    --rg-lb lib1 \
    --rg-sm sample1
```

### 跳过Fastp质控

```bash
python ngs_minimap2.py \
    -i sample_R1.fq.gz \
    -j sample_R2.fq.gz \
    -o output_dir \
    -r reference.fa \
    --skip-fastp \
    --rg-id sample1 \
    --rg-lb lib1 \
    --rg-sm sample1
```

### 仅保留唯一比对

```bash
python ngs_minimap2.py \
    -i sample_R1.fq.gz \
    -j sample_R2.fq.gz \
    -o output_dir \
    -r reference.fa \
    --uniq-mapping \
    --rg-id sample1 \
    --rg-lb lib1 \
    --rg-sm sample1
```

## 参数说明

### 必需参数

- `-i, --input1`: R1 fastq文件（或单端输入）
- `-o, --output_dir`: 输出目录
- `-r, --reference`: 参考基因组（FASTA格式，已索引）
- `--rg-id`: Read group ID
- `--rg-lb`: Read group文库
- `--rg-sm`: Read group样本名

### 可选参数

- `-j, --input2`: R2 fastq文件（双端模式）
- `-t, --thread`: 线程数（默认：8）
- `-p, --prefix`: 输出前缀（默认：sample）
- `--rg-pl`: Read group平台（默认：ILLUMINA）
- `-f, --output-format`: 输出格式，cram或bam（默认：cram）

### 标志选项

- `--dry-run`: 仅打印命令不执行
- `--keep-intermediate`: 保留中间文件
- `--no-markdup`: 跳过去重步骤
- `--uniq-mapping`: 过滤非唯一比对（samtools view -q 20 -F 3332）
- `--skip-fastp`: 跳过Fastp质控，直接使用原始fastq比对

## 输出文件

### Fastp步骤（如果未跳过）

- `{prefix}_fastp.html`: HTML格式的质控报告
- `{prefix}_fastp.json`: JSON格式的质控报告
- `{prefix}_clean_1.fq.gz`: 质控后的R1 reads
- `{prefix}_clean_2.fq.gz`: 质控后的R2 reads（仅双端）

### 比对步骤

- `{prefix}.sorted.bam`: 排序后的BAM文件
- `{prefix}.sorted.bam.csi`: BAM索引
- `{prefix}_minimap2.log`: Minimap2日志文件

### 去重步骤（如果启用）

- `{prefix}.rmdup.bam`: 去重后的BAM文件
- `{prefix}.rmdup.bam.bai`: BAM索引

### 最终输出

- CRAM模式：`{prefix}_sort.cram` 和 `{prefix}_sort.cram.crai`
- BAM模式：最终BAM文件及索引

## 示例

### 标准流程（包含所有步骤）

```bash
python ngs_minimap2.py \
    -i R1.fq.gz \
    -j R2.fq.gz \
    -o results \
    -r genome.fa \
    -t 32 \
    -p my_sample \
    --rg-id my_sample \
    --rg-lb my_library \
    --rg-sm my_sample \
    --output-format cram
```

### 快速比对（无质控、无去重、BAM输出）

```bash
python ngs_minimap2.py \
    -i R1.fq.gz \
    -j R2.fq.gz \
    -o results \
    -r genome.fa \
    --skip-fastp \
    --no-markdup \
    --output-format bam \
    --rg-id my_sample \
    --rg-lb my_library \
    --rg-sm my_sample
```

## 注意事项

1. 参考基因组必须已索引（`samtools faidx reference.fa`）
2. 去重功能需要安装Sambamba
3. 除非指定`--keep-intermediate`，否则会自动清理中间文件
4. 使用`--dry-run`可在执行前预览命令

## 更新日志

### v2.0

- 添加`--skip-fastp`选项，支持跳过质控步骤
- 添加`--uniq-mapping`选项，支持唯一比对过滤
- 改进中间文件清理逻辑
- 优化线程分配策略
- 完善错误处理和日志输出
