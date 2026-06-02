# TE注释流程快速使用指南

## 🚀 快速开始

### 1. 准备工作

```bash
# 确保基因组文件存在
ls genome.fasta

# 确保所有软件路径已配置（编辑TE_annotation.py顶部）
```

### 2. 运行流程

#### 方式A: 完整流程（推荐新手）

```bash
# 一键运行所有步骤
python TE_annotation.py -g genome.fasta -s MySpecies -t 16
```

#### 方式B: 分步运行（适合调试）

```bash
# Step 1: MITE-Hunter (约2-4小时)
python TE_annotation.py -g genome.fasta -s MySpecies --step1

# Step 2: LTR注释 (约4-8小时)
python TE_annotation.py -g genome.fasta -s MySpecies --step2

# Step 3: RepeatModeler + RepeatMasker (约8-12小时)
python TE_annotation.py -g genome.fasta -s MySpecies --step3 \
    --mite-lib 01_MITE/MySpecies.MITE.fa \
    --ltr-lib 04_LTR_Retriever/MySpecies.LTR.lib.fa
```

### 3. 查看结果

```bash
# 查看日志
tail -f MySpecies.TE_annotation.log

# 查看最终TE库
ls -lh 08_TE_Library/MySpecies.final.lib

# 查看RepeatMasker结果
head 07_RepeatMasker/genome.fasta.out
```

## 📊 运行时间参考

| 步骤                  | 小基因组 (≤1Gb) | 中等基因组 (1-3Gb) | 大基因组 (≥3Gb) |
| --------------------- | --------------- | ------------------ | --------------- |
| Step 1: MITE-Hunter   | 1-2小时         | 2-4小时            | 4-8小时         |
| Step 2: LTR注释       | 2-4小时         | 4-8小时            | 8-16小时        |
| Step 3: RepeatModeler | 4-8小时         | 8-16小时           | 16-32小时       |
| **总计**              | **7-14小时**    | **14-28小时**      | **28-56小时**   |

_注: 基于16线程，实际时间取决于基因组重复序列含量和硬件性能_

## 🎯 核心输出文件

| 文件                   | 用途                 | 下游分析               |
| ---------------------- | -------------------- | ---------------------- |
| `MySpecies.final.lib`  | **最终TE库**         | TE分类、比较基因组学   |
| `genome.fasta.out`     | RepeatMasker表格结果 | TE密度分析、基因组注释 |
| `genome.fasta.out.gff` | GFF格式TE注释        | 基因组浏览器可视化     |

## 🔧 常见问题快速解决

### 问题1: 找不到软件

```bash
# 检查软件路径
which ltr_finder
which RepeatMasker

# 编辑脚本，修改软件路径
vim TE_annotation.py +10  # 跳转到第10行（软件路径配置）
```

### 问题2: 内存不足

```bash
# 减少线程数
python TE_annotation.py -g genome.fasta -s MySpecies -t 8  # 从16减到8

# 或分步运行，每次只运行一个步骤
python TE_annotation.py -g genome.fasta -s MySpecies --step1
```

### 问题3: 中断后继续运行

```bash
# 脚本会自动跳过已完成的步骤
# 直接重新运行即可
python TE_annotation.py -g genome.fasta -s MySpecies -t 16
```

### 问题4: 检查中间结果

```bash
# 查看MITE库
grep -c "^>" 01_MITE/MySpecies.MITE.fa

# 查看LTR库
grep -c "^>" 04_LTR_Retriever/MySpecies.LTR.lib.fa

# 查看RepeatMasker覆盖度
awk '{sum+=$3-$2} END {print "Total TE length:", sum}' \
    07_RepeatMasker/genome.fasta.out
```

## 📝 参数说明

### 必需参数

| 参数            | 说明           | 示例              |
| --------------- | -------------- | ----------------- |
| `-g, --genome`  | 输入基因组文件 | `-g genome.fasta` |
| `-s, --species` | 物种名称       | `-s MySpecies`    |

### 可选参数

| 参数            | 说明         | 默认值 |
| --------------- | ------------ | ------ |
| `-t, --threads` | 线程数       | 8      |
| `--all`         | 运行所有步骤 | 默认   |
| `--step1`       | 只运行Step 1 | -      |
| `--step2`       | 只运行Step 2 | -      |
| `--step3`       | 只运行Step 3 | -      |
| `--mite-lib`    | MITE库文件   | -      |
| `--ltr-lib`     | LTR库文件    | -      |

## 🎓 学习资源

### 推荐阅读顺序

1. **基础概念**: `TE_annotation_README.md` → 概述章节
2. **流程详解**: `TE_annotation_README.md` → 流程图和步骤说明
3. **参数调优**: `TE_annotation_README.md` → 参数调优章节
4. **故障排除**: `TE_annotation_README.md` → 常见问题章节

### 参考文献

流程中使用的每个工具都有相应文献，详见 `TE_annotation_README.md` 的引用文献章节。

## 📧 获取帮助

```bash
# 查看帮助信息
python TE_annotation.py -h

# 查看示例
cat run_TE_annotation.sh
```

---

**祝您使用顺利！** 🎉
