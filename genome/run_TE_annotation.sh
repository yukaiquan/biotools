#!/bin/bash
# TE注释流程示例脚本

# ===================== 配置参数 =====================
GENOME="genome.fasta"          # 输入基因组文件（TR softmasked）
SPECIES="Species_name"         # 物种名称
THREADS=16                     # 线程数

# ===================== 运行流程 =====================

echo "=========================================="
echo "TE注释流程示例"
echo "=========================================="
echo "基因组: $GENOME"
echo "物种: $SPECIES"
echo "线程数: $THREADS"
echo "=========================================="

# 方式1: 运行完整流程（推荐）
echo "[方式1] 运行完整流程..."
python TE_annotation.py \
    -g $GENOME \
    -s $SPECIES \
    -t $THREADS \
    --all

# 方式2: 分步运行
# echo "[方式2] 分步运行..."
# 
# # Step 1: MITE-Hunter
# echo "[Step 1] MITE-Hunter..."
# python TE_annotation.py -g $GENOME -s $SPECIES --step1
# 
# # Step 2: LTR注释
# echo "[Step 2] LTR注释..."
# python TE_annotation.py -g $GENOME -s $SPECIES --step2
# 
# # Step 3: RepeatModeler/RepeatMasker
# echo "[Step 3] RepeatModeler/RepeatMasker..."
# python TE_annotation.py \
#     -g $GENOME \
#     -s $SPECIES \
#     --step3 \
#     --mite-lib 01_MITE/${SPECIES}.MITE.fa \
#     --ltr-lib 04_LTR_Retriever/${SPECIES}.LTR.lib.fa

echo "=========================================="
echo "TE注释流程完成！"
echo "=========================================="

# ===================== 结果检查 =====================
echo ""
echo "结果文件检查:"
echo ""

# 检查各步骤输出
check_file() {
    if [ -f "$1" ]; then
        size=$(du -h "$1" | cut -f1)
        echo "✅ $1 ($size)"
    else
        echo "❌ $1 (不存在)"
    fi
}

echo "Step 1输出:"
check_file "01_MITE/${SPECIES}.MITE.fa"

echo ""
echo "Step 2输出:"
check_file "04_LTR_Retriever/${SPECIES}.LTR.lib.fa"

echo ""
echo "Step 3输出:"
check_file "08_TE_Library/${SPECIES}.TE.lib"
check_file "08_TE_Library/${SPECIES}.final.lib"
check_file "07_RepeatMasker/${GENOME}.out"

echo ""
echo "日志文件:"
check_file "${SPECIES}.TE_annotation.log"

# ===================== 结果统计 =====================
echo ""
echo "=========================================="
echo "TE库统计:"
echo "=========================================="

# 统计各库序列数量
count_seqs() {
    if [ -f "$1" ]; then
        count=$(grep -c "^>" "$1" 2>/dev/null || echo "0")
        echo "$1: $count 条序列"
    fi
}

count_seqs "01_MITE/${SPECIES}.MITE.fa"
count_seqs "04_LTR_Retriever/${SPECIES}.LTR.lib.fa"
count_seqs "08_TE_Library/${SPECIES}.TE.lib"
count_seqs "08_TE_Library/${SPECIES}.final.lib"

echo ""
echo "完成时间: $(date)"
