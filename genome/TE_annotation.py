#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed 7 18 15:07:59 2023
Updated on 2026-05-13
TE注释流程（三步骤）：
1. MITE-Hunter搜索MITE转座子
2. LTR_Finder + LTR_Harvest + LTR_Retriever构建LTR库
3. RepeatModeler de novo搜索 + RepeatMasker全基因组注释

参考文献：
- MITE-Hunter: Han Y et al; 2010
- LTR_Finder: Xu Z at el; 2007
- LTR_Harvest: D. Ellinghaus et al; 2008
- LTR_Retriever: Ou S. at el; 2018
- RepeatModeler: Bedell, J.A. et al; 2000
- TEclass: György Abrusán at el; 2009
- Repbase: J, J. et al; 2005
- RepeatMasker: Smit, AFA et al; 2013
"""
import os
import sys
import argparse
import logging
import subprocess
import shutil
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import List, Tuple, Optional
from datetime import datetime

# ===================== 软件路径配置 =====================
# MITE-Hunter
MITEHUNTER = "/publicssd/share/h13713/soft/MITEHunter/MITE_Hunter_manager.pl"

# LTR相关工具
LTR_FINDER = "/publicssd/share/h13713/soft/ltr_finder"
LTR_HARVEST = "/public/software/apps/geta/biosoft/genometools-1.6.2/bin/gt"  # GenomeTools
LTR_RETRIEVER = "/publicssd/share/h13713/soft/LTR_retriever/LTR_retriever"

# RepeatModeler
REPEATMODELER = "/public/software/apps/geta/biosoft/RepeatModeler-2.0.5/RepeatModeler"
BUILD_DATABASE = "/public/software/apps/geta/biosoft/RepeatModeler-2.0.5/BuildDatabase"

# TEclass
TESORTED = "/publicssd/share/h13713/soft"
TECLASS = "/public/share/acfaa2ssz7/soft/TEclass/TEclassTest.pl"

# RepeatMasker
REPEATMASKER = "/public/software/apps/geta/biosoft/RepeatMasker-4.1.6/RepeatMasker"
REPEATMASKER_LIB = "/public/share/acfaa2ssz7/soft/RepBase/RepBase.lib"

# 其他工具
SAMTOOLS = "samtools"
SEQTK = "seqtk"
CD_HIT = "cd-hit-est"

# ===================== 输出目录配置 =====================
MITEHUNTER_OUTDIR = "03_MITE"
LTR_FINDER_OUTDIR = "04_LTR_Finder"
LTR_HARVEST_OUTDIR = "05_LTR_Harvest"
LTR_RETRIEVER_OUTDIR = "06_LTR_Retriever"
REPEATMODELER_OUTDIR = "07_RepeatModeler"
TE_CLASS_OUTDIR = "08_TEclass"
REPEATMASKER_OUTDIR = "09_RepeatMasker"
TE_LIBRARY_OUTDIR = "10_TE_Library"

# ===================== 全局环境变量设置 =====================
tool_paths = [
    os.path.dirname(MITEHUNTER),
    os.path.dirname(LTR_FINDER),
    os.path.dirname(LTR_RETRIEVER),
    os.path.dirname(REPEATMODELER),
    os.path.dirname(TECLASS),
    os.path.dirname(REPEATMASKER),
]

for tool_path in tool_paths:
    if tool_path and os.path.exists(tool_path):
        os.environ["PATH"] = f"{tool_path}:{os.environ.get('PATH', '')}"

CWD = os.getcwd()

# ===================== 日志配置 =====================
def setup_logging(log_file: str = "TE_annotation.log") -> None:
    """配置日志记录"""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )


def get_timestamp() -> str:
    """获取当前时间戳"""
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


# ===================== 辅助函数 =====================
def run_cmd(cmd: str, description: str = "", check: bool = True) -> int:
    """
    运行命令并等待完成
    
    Args:
        cmd: 要执行的命令
        description: 命令描述
        check: 是否检查返回码
    
    Returns:
        返回码
    """
    if cmd == "NA":
        logging.info(f"跳过步骤: {description}")
        return 0
    
    logging.info(f"执行: {description}")
    logging.info(f"命令: {cmd}")
    print(f"\n[{get_timestamp()}] 执行: {description}")
    print(f"命令: {cmd}\n")
    
    try:
        result = subprocess.run(
            cmd,
            shell=True,
            check=check,
            capture_output=True,
            text=True,
            env=os.environ.copy()
        )
        
        if result.stdout:
            logging.info(f"标准输出:\n{result.stdout}")
        if result.stderr:
            # 某些工具会将正常输出输出到stderr
            if "error" in result.stderr.lower() or "fail" in result.stderr.lower():
                logging.warning(f"标准错误:\n{result.stderr}")
            else:
                logging.info(f"标准错误:\n{result.stderr}")
        
        return result.returncode
    except subprocess.CalledProcessError as e:
        logging.error(f"命令执行失败！返回码: {e.returncode}")
        logging.error(f"错误输出:\n{e.stderr}")
        if check:
            sys.exit(1)
        return e.returncode


def check_file_exists(file_path: str, description: str = "") -> bool:
    """检查文件是否存在"""
    if os.path.exists(file_path):
        logging.info(f"{description} 文件已存在: {file_path}")
        return True
    return False


def mkdir(directory: str) -> None:
    """创建目录"""
    if not os.path.exists(directory):
        os.makedirs(directory)
        logging.info(f"创建目录: {directory}")


def get_genome_size(genome: str) -> int:
    """获取基因组大小"""
    fai_file = f"{genome}.fai"
    if not os.path.exists(fai_file):
        run_cmd(f"samtools faidx {genome}", "创建基因组索引")
    
    genome_size = 0
    with open(fai_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            genome_size += int(parts[1])
    return genome_size


def get_chromosomes(genome: str) -> List[str]:
    """从基因组索引文件获取染色体名称"""
    fai_file = f"{genome}.fai"
    if not os.path.exists(fai_file):
        run_cmd(f"samtools faidx {genome}", "创建基因组索引")
    
    chroms = []
    with open(fai_file, 'r') as f:
        for line in f:
            chrom = line.strip().split('\t')[0]
            chroms.append(chrom)
    return chroms


# ===================== Step 1: MITE-Hunter =====================
class MITEHunter:
    """MITE-Hunter搜索MITE转座子"""
    
    def __init__(self, genome: str, species: str, threads: int = 8):
        self.genome = os.path.abspath(genome)
        self.species = species
        self.threads = threads
        self.outdir = MITEHUNTER_OUTDIR
        self.output_fasta = os.path.join(self.outdir, f"{species}.MITE.fa")
    
    def calculate_sample_ratio(self) -> float:
        """
        计算MITE-Hunter的采样比例
        对于700Mb以下的基因组用1.0
        参数可以设置为 1/(实际基因组大小/700)
        """
        genome_size_mb = get_genome_size(self.genome) / 1_000_000
        ratio = 1.0 / (genome_size_mb / 700.0)
        ratio = max(0.01, min(1.0, ratio))  # 限制在0.01-1.0之间
        logging.info(f"基因组大小: {genome_size_mb:.2f} Mb, MITE-Hunter采样比例: {ratio:.4f}")
        return ratio
    
    def run(self) -> bool:
        """运行MITE-Hunter"""
        mkdir(self.outdir)
        os.chdir(self.outdir)
        
        # 检查输出文件是否存在
        if check_file_exists(self.output_fasta, "MITE-Hunter输出"):
            return True
        
        # 计算采样比例
        sample_ratio = self.calculate_sample_ratio()
        
        # 构建命令
        # perl MITE_Hunter_manager.pl -i input.fa -S 12345678 -n 20 -P 0.22 -c 16 -g species_name
        cmd = (
            f"perl {MITEHUNTER} "
            f"-i {self.genome} "
            f"-S 12345678 "
            f"-n 20 "
            f"-P {sample_ratio:.4f} "
            f"-c {self.threads} "
            f"-g {self.species}"
        )
        
        returncode = run_cmd(cmd, "MITE-Hunter搜索MITE转座子", check=False)
        
        # 检查结果文件
        if os.path.exists(self.output_fasta):
            logging.info(f"MITE-Hunter完成，输出文件: {self.output_fasta}")
            return True
        else:
            logging.error("MITE-Hunter未生成输出文件")
            return False
    
    def get_output(self) -> str:
        """获取MITE-Hunter输出文件路径"""
        return self.output_fasta


# ===================== Step 2: LTR注释 =====================
class LTRAnnotation:
    """LTR_Finder + LTR_Harvest + LTR_Retriever构建LTR库"""
    
    def __init__(self, genome: str, species: str, threads: int = 8):
        self.genome = os.path.abspath(genome)
        self.species = species
        self.threads = threads
        self.ltr_finder_outdir = LTR_FINDER_OUTDIR
        self.ltr_harvest_outdir = LTR_HARVEST_OUTDIR
        self.ltr_retriever_outdir = LTR_RETRIEVER_OUTDIR
        self.output_fasta = os.path.join(self.ltr_retriever_outdir, f"{species}.LTR.lib.fa")
    
    def run_ltr_finder_parallel(self) -> bool:
        """并行运行LTR_Finder（按染色体）"""
        mkdir(self.ltr_finder_outdir)
        os.chdir(self.ltr_finder_outdir)
        
        # 获取染色体列表
        chroms = get_chromosomes(self.genome)
        logging.info(f"检测到 {len(chroms)} 条染色体")
        
        # 检查是否所有染色体都已完成
        all_scn_files = [f"{self.species}.ltr.finder_{chrom}.scn" for chrom in chroms]
        if all(check_file_exists(f, f"LTR_Finder结果 ({chrom})") for f in all_scn_files):
            logging.info("所有染色体的LTR_Finder结果已存在，跳过")
            return True
        
        # 构建命令列表
        cmd_list = []
        for chrom in chroms:
            output_fasta = f"{self.species}.ssrMasked_{chrom}.fasta"
            output_scn = f"{self.species}.ltr.finder_{chrom}.scn"
            
            # 提取染色体序列 + LTR_Finder
            cmd = (
                f"samtools faidx {self.genome} {chrom} > {output_fasta} && "
                f"{LTR_FINDER} -D 15000 -d 1000 -L 7000 -l 100 -p {self.threads} -C -M 0.9 "
                f"{output_fasta} > {output_scn}"
            )
            cmd_list.append((chrom, cmd))
        
        # 并行执行
        logging.info(f"使用 {self.threads} 个线程并行运行LTR_Finder...")
        with ProcessPoolExecutor(max_workers=self.threads) as executor:
            futures = {
                executor.submit(run_cmd, cmd, f"LTR_Finder ({chrom})", check=False): chrom
                for chrom, cmd in cmd_list
            }
            
            completed = 0
            for future in as_completed(futures):
                chrom = futures[future]
                try:
                    future.result()
                    completed += 1
                    logging.info(f"LTR_Finder完成: {chrom} ({completed}/{len(chroms)})")
                except Exception as e:
                    logging.error(f"LTR_Finder失败: {chrom}, 错误: {str(e)}")
        
        return True
    
    def run_ltr_harvest(self) -> bool:
        """运行LTR_Harvest"""
        mkdir(self.ltr_harvest_outdir)
        os.chdir(self.ltr_harvest_outdir)
        
        output_scn = f"{self.species}.harvest.scn"
        
        # 检查输出文件是否存在
        if check_file_exists(output_scn, "LTR_Harvest输出"):
            return True
        
        # 创建基因组软链接
        genome_basename = os.path.basename(self.genome)
        if not os.path.exists(genome_basename):
            os.symlink(self.genome, genome_basename)
        
        # Step 1: gt suffixerator - 构建索引
        cmd_suffix = (
            f"gt suffixerator -db {genome_basename} "
            f"-indexname {self.species} "
            f"-tis -suf -lcp -des -ssp -sds -dna"
        )
        run_cmd(cmd_suffix, "LTR_Harvest构建索引")
        
        # Step 2: gt ltrharvest - 搜索LTR
        cmd_harvest = (
            f"gt ltrharvest -index {self.species} "
            f"-similar 90 -vic 10 -seed 20 -seqids yes "
            f"-minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 "
            f"-motif TGCA -motifmis 1 > {output_scn}"
        )
        run_cmd(cmd_harvest, "LTR_Harvest搜索LTR")
        
        return True
    
    def merge_ltr_finder_results(self) -> str:
        """合并LTR_Finder结果"""
        os.chdir(self.ltr_finder_outdir)
        
        merged_scn = f"{self.species}.ltr.finder.scn"
        
        # 检查合并文件是否存在
        if check_file_exists(merged_scn, "LTR_Finder合并结果"):
            return os.path.join(self.ltr_finder_outdir, merged_scn)
        
        # 获取所有.scn文件
        chroms = get_chromosomes(self.genome)
        scn_files = [f"{self.species}.ltr.finder_{chrom}.scn" for chrom in chroms]
        
        # 合并
        with open(merged_scn, 'w') as out_f:
            for scn_file in scn_files:
                if os.path.exists(scn_file):
                    with open(scn_file, 'r') as in_f:
                        out_f.write(in_f.read())
                else:
                    logging.warning(f"LTR_Finder结果文件不存在: {scn_file}")
        
        logging.info(f"合并LTR_Finder结果: {merged_scn}")
        return os.path.join(self.ltr_finder_outdir, merged_scn)
    
    def run_ltr_retriever(self) -> bool:
        """运行LTR_Retriever构建LTR库"""
        mkdir(self.ltr_retriever_outdir)
        os.chdir(self.ltr_retriever_outdir)
        
        # 检查输出文件是否存在
        if check_file_exists(self.output_fasta, "LTR_Retriever输出"):
            return True
        
        # 获取LTR_Finder和LTR_Harvest的结果
        ltr_finder_scn = self.merge_ltr_finder_results()
        ltr_harvest_scn = os.path.join(CWD, self.ltr_harvest_outdir, f"{self.species}.harvest.scn")
        
        # LTR_Retriever命令
        cmd = (
            f"perl {LTR_RETRIEVER} "
            f"-genome {self.genome} "
            f"-ltrfinder {ltr_finder_scn} "
            f"-ltrharvest {ltr_harvest_scn} "
            f"-threads {self.threads}"
        )
        
        run_cmd(cmd, "LTR_Retriever构建LTR库", check=False)
        
        # 检查结果
        if os.path.exists(self.output_fasta):
            logging.info(f"LTR_Retriever完成，输出文件: {self.output_fasta}")
            return True
        else:
            logging.error("LTR_Retriever未生成输出文件")
            return False
    
    def run(self) -> bool:
        """运行完整的LTR注释流程"""
        logging.info("=" * 80)
        logging.info("Step 2: LTR注释流程")
        logging.info("=" * 80)
        
        # Step 2.1: LTR_Finder
        logging.info("\n[2.1] 运行LTR_Finder...")
        if not self.run_ltr_finder_parallel():
            return False
        
        # Step 2.2: LTR_Harvest
        logging.info("\n[2.2] 运行LTR_Harvest...")
        if not self.run_ltr_harvest():
            return False
        
        # Step 2.3: LTR_Retriever
        logging.info("\n[2.3] 运行LTR_Retriever...")
        if not self.run_ltr_retriever():
            return False
        
        logging.info("\nStep 2: LTR注释流程完成！")
        return True
    
    def get_output(self) -> str:
        """获取LTR库输出文件路径"""
        return self.output_fasta


# ===================== Step 3: RepeatModeler + RepeatMasker =====================
class RepeatAnnotation:
    """RepeatModeler de novo搜索 + RepeatMasker全基因组注释"""
    
    def __init__(self, genome: str, species: str, threads: int = 8):
        self.genome = os.path.abspath(genome)
        self.species = species
        self.threads = threads
        self.te_lib_dir = TE_LIBRARY_OUTDIR
        self.repeatmodeler_outdir = REPEATMODELER_OUTDIR
        self.teclass_outdir = TE_CLASS_OUTDIR
        self.repeatmasker_outdir = REPEATMASKER_OUTDIR
        
        # 输出文件
        self.te_lib = os.path.join(self.te_lib_dir, f"{species}.TE.lib")
        self.repmod_lib = os.path.join(self.repeatmodeler_outdir, f"{species}.RepMod.lib")
        self.final_lib = os.path.join(self.te_lib_dir, f"{species}.final.lib")
    
    def build_te_library(self, mite_lib: str, ltr_lib: str) -> str:
        """构建TE库（整合MITE和LTR）"""
        mkdir(self.te_lib_dir)
        os.chdir(self.te_lib_dir)
        
        if check_file_exists(self.te_lib, "TE库"):
            return self.te_lib
        
        logging.info("构建TE库（MITE + LTR）...")
        
        # 合并MITE和LTR库
        with open(self.te_lib, 'w') as out_f:
            # 写入MITE序列
            if os.path.exists(mite_lib):
                logging.info(f"添加MITE库: {mite_lib}")
                with open(mite_lib, 'r') as in_f:
                    out_f.write(in_f.read())
            
            # 写入LTR序列
            if os.path.exists(ltr_lib):
                logging.info(f"添加LTR库: {ltr_lib}")
                with open(ltr_lib, 'r') as in_f:
                    out_f.write(in_f.read())
        
        # 去重（可选，使用CD-HIT）
        # cd-hit-est -i TE.lib -o TE.lib.nr -c 0.9 -n 8
        
        logging.info(f"TE库构建完成: {self.te_lib}")
        return self.te_lib
    
    def hardmask_genome(self, genome: str, te_lib: str) -> str:
        """使用TE库对基因组进行hardmask"""
        mkdir(self.repeatmasker_outdir)
        os.chdir(self.repeatmasker_outdir)
        
        output_genome = f"{self.species}.hardmasked.fasta"
        
        if check_file_exists(output_genome, "Hardmask基因组"):
            return os.path.join(self.repeatmasker_outdir, output_genome)
        
        logging.info("使用TE库对基因组进行hardmask...")
        
        # RepeatMasker hardmask
        cmd = (
            f"{REPEATMASKER} -pa {self.threads} "
            f"-lib {te_lib} "
            f"-xsmall "
            f"-dir . "
            f"{genome}"
        )
        
        run_cmd(cmd, "RepeatMasker hardmask")
        
        # RepeatMasker输出为 .masked 文件
        masked_file = genome.replace(".fasta", ".fasta.masked")
        if os.path.exists(masked_file):
            shutil.move(masked_file, output_genome)
        
        logging.info(f"Hardmask基因组: {output_genome}")
        return os.path.join(self.repeatmasker_outdir, output_genome)
    
    def run_repeatmodeler(self, genome: str) -> str:
        """运行RepeatModeler进行de novo搜索"""
        mkdir(self.repeatmodeler_outdir)
        os.chdir(self.repeatmodeler_outdir)
        
        if check_file_exists(self.repmod_lib, "RepeatModeler输出"):
            return self.repmod_lib
        
        logging.info("运行RepeatModeler进行de novo重复序列搜索...")
        
        # Step 1: BuildDatabase
        cmd_build = (
            f"{BUILD_DATABASE} -name {self.species} "
            f"-engine ncbi {genome}"
        )
        run_cmd(cmd_build, "RepeatModeler构建数据库")
        
        # Step 2: RepeatModeler
        cmd_repeatmodeler = (
            f"{REPEATMODELER} -database {self.species} "
            f"-pa {self.threads} "
            f"-LTRStruct"
        )
        run_cmd(cmd_repeatmodeler, "RepeatModeler de novo搜索")
        
        # RepeatModeler输出文件通常为 families-classified.stk 或 consensi.fa.classified
        # 需要找到并复制
        possible_outputs = [
            "families-classified.stk",
            "consensi.fa.classified",
            f"RM_{self.species}.consensi.fa.classified"
        ]
        
        repmod_output = None
        for output in possible_outputs:
            if os.path.exists(output):
                repmod_output = output
                break
        
        if repmod_output:
            # 转换为FASTA格式（如果是stk格式）
            if repmod_output.endswith('.stk'):
                # 使用Alignment格式转换
                fasta_output = repmod_output.replace('.stk', '.fa')
                cmd_convert = f"perl -ne 'if(/^>(\\S+)/){{print \">$1\\n\"}} elsif(/^([^#\\n])/){{print \"$1\\n\"}}' {repmod_output} > {fasta_output}"
                run_cmd(cmd_convert, "转换stk为FASTA格式")
                shutil.copy(fasta_output, self.repmod_lib)
            else:
                shutil.copy(repmod_output, self.repmod_lib)
            
            logging.info(f"RepeatModeler输出: {self.repmod_lib}")
        else:
            logging.warning("未找到RepeatModeler输出文件")
            return ""
        
        return self.repmod_lib
    
    def classify_unknown_te(self, te_lib: str) -> str:
        """使用TEclass对Unknown重复序列进行分类"""
        mkdir(self.teclass_outdir)
        os.chdir(self.teclass_outdir)
        
        classified_lib = f"{self.species}.classified.lib"
        
        if check_file_exists(classified_lib, "TEclass分类结果"):
            return os.path.join(self.teclass_outdir, classified_lib)
        
        logging.info("使用TEclass对Unknown重复序列进行分类...")
        
        # TEclass命令
        cmd = f"perl {TECLASS} -i {te_lib} -o {classified_lib}"
        run_cmd(cmd, "TEclass分类")
        
        return os.path.join(self.teclass_outdir, classified_lib)
    
    def build_final_library(self, te_lib: str, repmod_lib: str, repbase_lib: str = None) -> str:
        """构建最终库文件（整合所有库）"""
        mkdir(self.te_lib_dir)
        os.chdir(self.te_lib_dir)
        
        if check_file_exists(self.final_lib, "最终TE库"):
            return self.final_lib
        
        logging.info("构建最终TE库（TE.lib + RepMod.lib + Repbase）...")
        
        with open(self.final_lib, 'w') as out_f:
            # 写入TE.lib
            if os.path.exists(te_lib):
                logging.info(f"添加TE库: {te_lib}")
                with open(te_lib, 'r') as in_f:
                    out_f.write(in_f.read())
            
            # 写入RepMod.lib
            if os.path.exists(repmod_lib):
                logging.info(f"添加RepMod库: {repmod_lib}")
                with open(repmod_lib, 'r') as in_f:
                    out_f.write(in_f.read())
            
            # 写入Repbase（如果存在）
            if repbase_lib and os.path.exists(repbase_lib):
                logging.info(f"添加Repbase库: {repbase_lib}")
                with open(repbase_lib, 'r') as in_f:
                    out_f.write(in_f.read())
        
        # 去重
        nr_lib = self.final_lib + ".nr"
        cmd_nr = f"{CD_HIT} -i {self.final_lib} -o {nr_lib} -c 0.9 -n 8"
        run_cmd(cmd_nr, "CD-HIT去重")
        
        # 使用去重后的库
        if os.path.exists(nr_lib):
            shutil.move(nr_lib, self.final_lib)
        
        logging.info(f"最终TE库: {self.final_lib}")
        return self.final_lib
    
    def run_repeatmasker(self, genome: str, te_lib: str) -> bool:
        """使用RepeatMasker对全基因组进行重复序列搜索"""
        mkdir(self.repeatmasker_outdir)
        os.chdir(self.repeatmasker_outdir)
        
        logging.info("使用最终TE库对全基因组进行RepeatMasker注释...")
        
        # RepeatMasker命令
        cmd = (
            f"{REPEATMASKER} -pa {self.threads} "
            f"-lib {te_lib} "
            f"-species {self.species} "
            f"-xsmall "
            f"-dir . "
            f"-gff "
            f"{genome}"
        )
        
        run_cmd(cmd, "RepeatMasker全基因组注释")
        
        logging.info("RepeatMasker完成！")
        return True
    
    def run(self, mite_lib: str, ltr_lib: str) -> bool:
        """运行完整的RepeatModeler + RepeatMasker流程"""
        logging.info("=" * 80)
        logging.info("Step 3: RepeatModeler + RepeatMasker注释")
        logging.info("=" * 80)
        
        # Step 3.1: 构建TE库
        logging.info("\n[3.1] 构建TE库（MITE + LTR）...")
        te_lib = self.build_te_library(mite_lib, ltr_lib)
        
        # Step 3.2: Hardmask基因组
        logging.info("\n[3.2] 对基因组进行hardmask...")
        hardmasked_genome = self.hardmask_genome(self.genome, te_lib)
        
        # Step 3.3: RepeatModeler
        logging.info("\n[3.3] 运行RepeatModeler...")
        repmod_lib = self.run_repeatmodeler(hardmasked_genome)
        
        # Step 3.4: TEclass分类
        logging.info("\n[3.4] 使用TEclass分类Unknown序列...")
        classified_lib = self.classify_unknown_te(repmod_lib)
        
        # Step 3.5: 构建最终库
        logging.info("\n[3.5] 构建最终TE库...")
        final_lib = self.build_final_library(te_lib, classified_lib, REPEATMASKER_LIB)
        
        # Step 3.6: RepeatMasker全基因组注释
        logging.info("\n[3.6] RepeatMasker全基因组注释...")
        if not self.run_repeatmasker(self.genome, final_lib):
            return False
        
        logging.info("\nStep 3: RepeatModeler + RepeatMasker注释完成！")
        return True
    
    def get_final_library(self) -> str:
        """获取最终TE库文件路径"""
        return self.final_lib


# ===================== 主程序 =====================
def argparse_init():
    """初始化命令行参数解析器"""
    parser = argparse.ArgumentParser(
        description="TE注释流程（三步骤）：MITE-Hunter + LTR注释 + RepeatModeler/RepeatMasker",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  # 运行完整流程
  python TE_annotation.py -g genome.fasta -s species_name -t 16
  
  # 只运行Step 1 (MITE-Hunter)
  python TE_annotation.py -g genome.fasta -s species_name --step1
  
  # 只运行Step 2 (LTR注释)
  python TE_annotation.py -g genome.fasta -s species_name --step2
  
  # 只运行Step 3 (RepeatModeler/RepeatMasker)
  python TE_annotation.py -g genome.fasta -s species_name --step3 --mite-lib MITE.fa --ltr-lib LTR.fa
        """
    )
    
    parser.add_argument("-g", "--genome", required=True,
                        help="输入基因组FASTA文件（TR softmasked）")
    parser.add_argument("-s", "--species", required=True,
                        help="物种名称（用于输出文件命名）")
    parser.add_argument("-t", "--threads", type=int, default=8,
                        help="并行线程数（默认: 8）")
    
    # 步骤控制
    step_group = parser.add_mutually_exclusive_group()
    step_group.add_argument("--all", action="store_true",
                            help="运行所有步骤（默认）")
    step_group.add_argument("--step1", action="store_true",
                            help="只运行Step 1: MITE-Hunter")
    step_group.add_argument("--step2", action="store_true",
                            help="只运行Step 2: LTR注释")
    step_group.add_argument("--step3", action="store_true",
                            help="只运行Step 3: RepeatModeler/RepeatMasker")
    
    # Step 3需要的输入文件
    parser.add_argument("--mite-lib", 
                        help="MITE库文件（用于Step 3）")
    parser.add_argument("--ltr-lib",
                        help="LTR库文件（用于Step 3）")
    
    return parser.parse_args()


def main():
    """主程序"""
    args = argparse_init()
    
    # 设置日志
    log_file = f"{args.species}.TE_annotation.log"
    setup_logging(log_file)
    
    logging.info("=" * 80)
    logging.info("TE注释流程启动")
    logging.info(f"基因组文件: {args.genome}")
    logging.info(f"物种名称: {args.species}")
    logging.info(f"线程数: {args.threads}")
    logging.info("=" * 80)
    
    # 检查基因组文件
    if not os.path.exists(args.genome):
        logging.error(f"基因组文件不存在: {args.genome}")
        sys.exit(1)
    
    # 确定要运行的步骤
    run_step1 = args.step1 or args.all or (not args.step1 and not args.step2 and not args.step3)
    run_step2 = args.step2 or args.all or (not args.step1 and not args.step2 and not args.step3)
    run_step3 = args.step3 or args.all or (not args.step1 and not args.step2 and not args.step3)
    
    mite_lib = args.mite_lib
    ltr_lib = args.ltr_lib
    
    try:
        # Step 1: MITE-Hunter
        if run_step1:
            logging.info("\n" + "=" * 80)
            logging.info("Step 1: MITE-Hunter搜索MITE转座子")
            logging.info("=" * 80)
            
            mite = MITEHunter(args.genome, args.species, args.threads)
            if mite.run():
                mite_lib = mite.get_output()
                logging.info(f"Step 1完成！MITE库: {mite_lib}")
            else:
                logging.error("Step 1失败！")
                sys.exit(1)
        
        # Step 2: LTR注释
        if run_step2:
            logging.info("\n" + "=" * 80)
            logging.info("Step 2: LTR注释（LTR_Finder + LTR_Harvest + LTR_Retriever）")
            logging.info("=" * 80)
            
            ltr = LTRAnnotation(args.genome, args.species, args.threads)
            if ltr.run():
                ltr_lib = ltr.get_output()
                logging.info(f"Step 2完成！LTR库: {ltr_lib}")
            else:
                logging.error("Step 2失败！")
                sys.exit(1)
        
        # Step 3: RepeatModeler + RepeatMasker
        if run_step3:
            # 检查输入库文件
            if not mite_lib or not os.path.exists(mite_lib):
                logging.error("Step 3需要MITE库文件，请使用 --mite-lib 参数或先运行Step 1")
                sys.exit(1)
            
            if not ltr_lib or not os.path.exists(ltr_lib):
                logging.error("Step 3需要LTR库文件，请使用 --ltr-lib 参数或先运行Step 2")
                sys.exit(1)
            
            logging.info("\n" + "=" * 80)
            logging.info("Step 3: RepeatModeler + RepeatMasker注释")
            logging.info("=" * 80)
            
            repeat = RepeatAnnotation(args.genome, args.species, args.threads)
            if repeat.run(mite_lib, ltr_lib):
                final_lib = repeat.get_final_library()
                logging.info(f"Step 3完成！最终TE库: {final_lib}")
            else:
                logging.error("Step 3失败！")
                sys.exit(1)
        
        logging.info("\n" + "=" * 80)
        logging.info("TE注释流程全部完成！")
        logging.info("=" * 80)
        
        # 输出结果文件列表
        logging.info("\n结果文件:")
        if run_step1:
            logging.info(f"  MITE库: {mite_lib}")
        if run_step2:
            logging.info(f"  LTR库: {ltr_lib}")
        if run_step3:
            logging.info(f"  最终TE库: {final_lib}")
            logging.info(f"  RepeatMasker结果目录: {REPEATMASKER_OUTDIR}/")
        
        logging.info(f"\n日志文件: {log_file}")
        
    except KeyboardInterrupt:
        logging.error("\n用户中断执行！")
        sys.exit(1)
    except Exception as e:
        logging.error(f"\n执行出错: {str(e)}")
        import traceback
        logging.error(traceback.format_exc())
        sys.exit(1)


if __name__ == "__main__":
    main()
