#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed 7 18 15:07:59 2023
This script is used to annotate the genome structure of the species of interest.
"""
import os
import sys
import argparse
import logging
import subprocess


# ===================== 软件路径配置 =====================
GMATA = "/publicssd/share/h13713/soft/GMATA/gmata.pl"
GMATA_PATH = "/publicssd/share/h13713/soft/GMATA"
GMATA_CONFIG = "/publicssd/share/h13713/soft/GMATA/default_cfg.txt"
TRF = "/public/software/apps/geta/biosoft/RepeatMasker-4.1.6/trf"
GMATA_MASK = "/publicssd/share/h13713/soft/GMATA/gssrmsk.pl"
DTA2BED = "/publicssd/share/h13713/soft/dta2bed.py"

# ===================== 输出目录配置 =====================
GMATA_OUTDIR = "01_GMATA"
TRF_OUTDIR = "02_TRF"

# ===================== 全局环境变量设置（关键修改） =====================
# 将GMATA路径添加到系统PATH（让shell能找到GMATA的内部工具）
os.environ["PATH"] = f"{GMATA_PATH}:{os.environ.get('PATH', '')}"
# 将GMATA路径添加到Perl模块搜索路径（解决"Can't locate XXX.pm"错误）
os.environ["PERL5LIB"] = f"{GMATA_PATH}:{os.environ.get('PERL5LIB', '')}"

CWD = os.getcwd()


def run_cmd_wait(cmd: str) -> int:
    """运行命令并等待完成，检查返回码"""
    if cmd == "NA":
        logging.info("跳过已完成的步骤")
        return 0
    
    logging.info(f"执行命令: {cmd}")
    print(f"\n[执行] {cmd}")
    
    # 使用check=True让命令失败时自动抛出异常
    try:
        result = subprocess.run(
            cmd, 
            shell=True, 
            check=True,
            capture_output=True,
            text=True
        )
        if result.stdout:
            logging.info(f"标准输出:\n{result.stdout}")
        if result.stderr:
            logging.info(f"标准错误:\n{result.stderr}")
        return result.returncode
    except subprocess.CalledProcessError as e:
        logging.error(f"命令执行失败！返回码: {e.returncode}")
        logging.error(f"错误输出:\n{e.stderr}")
        print(f"\n[错误] 命令执行失败: {cmd}")
        print(f"错误信息: {e.stderr}")
        sys.exit(1)


class File:
    """文件操作类"""
    def __init__(self, infile: str):
        self.infile = infile

    def check_file(self) -> bool:
        """检查文件是否存在"""
        if not os.path.exists(self.infile):
            logging.error(f"文件不存在: {self.infile}")
            return False
        return True

    def mv_file(self, outdir: str) -> str:
        """移动文件到指定目录"""
        if not self.check_file():
            return "NA"
        dest_path = os.path.join(outdir, self.infile)
        if os.path.exists(dest_path):
            logging.info(f"目标文件已存在，跳过移动: {dest_path}")
            return "NA"
        return f"mv {self.infile} {outdir}"

    def mkdir(self) -> str:
        """创建目录"""
        if not os.path.exists(self.infile):
            return f"mkdir -p {self.infile}"
        return "NA"


class TR_Annotation:
    """串联重复序列注释类"""
    def __init__(self, genome: str, species: str):
        # 每次初始化回到运行目录
        os.chdir(CWD)
        self.genome = genome
        self.genome_basename = os.path.basename(genome)
        self.outdir = GMATA_OUTDIR
        self.trf_outdir = TRF_OUTDIR
        self.species = species
        
        # 输出文件路径
        self.gmata_out_file = f"{self.genome_basename}.ssr"
        self.mask_out_file = f"{self.genome_basename}.ssrMasked.fasta"
        self.trf_out_file = f"{self.mask_out_file}.2.7.7.80.10.50.500.dat"

    def run_gmata(self) -> str:
        """运行GMATA进行SSR预测"""
        os.chdir(os.path.join(CWD, self.outdir))
        output_path = os.path.join(CWD, self.outdir, self.gmata_out_file)
        
        if os.path.exists(output_path):
            logging.info(f"GMATA结果已存在，跳过: {output_path}")
            return "NA"
        
        # 使用绝对路径引用基因组文件，避免相对路径问题
        genome_abs_path = os.path.abspath(os.path.join(CWD, self.genome))
        return f"perl {GMATA} -c {GMATA_CONFIG} -i {genome_abs_path}"

    def gmata_mask(self) -> str:
        """使用GMATA结果屏蔽基因组中的SSR序列"""
        os.chdir(os.path.join(CWD, self.outdir))
        output_path = os.path.join(CWD, self.outdir, self.mask_out_file)
        
        if os.path.exists(output_path):
            logging.info(f"SSR屏蔽结果已存在，跳过: {output_path}")
            return "NA"
        
        genome_abs_path = os.path.abspath(os.path.join(CWD, self.genome))
        ssr_abs_path = os.path.abspath(os.path.join(CWD, self.outdir, self.gmata_out_file))
        return f"perl {GMATA_MASK} -i {genome_abs_path} -s 1 -r {ssr_abs_path}"

    def run_trf(self) -> str:
        """运行TRF进行串联重复序列预测"""
        os.chdir(os.path.join(CWD, self.trf_outdir))
        output_path = os.path.join(CWD, self.trf_outdir, self.trf_out_file)
        
        if os.path.exists(output_path):
            logging.info(f"TRF结果已存在，跳过: {output_path}")
            return "NA"
        
        mask_fasta_abs_path = os.path.abspath(os.path.join(CWD, self.outdir, self.mask_out_file))
        return f"{TRF} {mask_fasta_abs_path} 2 7 7 80 10 50 500 -d -h"

    def run_dta2bed(self) -> str:
        """将TRF的dat格式转换为bed格式"""
        os.chdir(os.path.join(CWD, self.trf_outdir))
        output_path = os.path.join(CWD, self.trf_outdir, f"{self.trf_out_file}.bed")
        
        if os.path.exists(output_path):
            logging.info(f"BED文件已存在，跳过: {output_path}")
            return "NA"
        
        return f"python3 {DTA2BED} --dat {self.trf_out_file} --bed {self.trf_out_file}.bed"

    def mask_2_lower(self) -> str:
        """使用seqtk将bed区域转换为小写字母"""
        os.chdir(os.path.join(CWD, self.trf_outdir))
        output_path = os.path.join(CWD, self.trf_outdir, f"{self.species}.ssrMasked.fasta")
        
        if os.path.exists(output_path):
            logging.info(f"最终屏蔽基因组已存在，跳过: {output_path}")
            return "NA"
        
        mask_fasta_abs_path = os.path.abspath(os.path.join(CWD, self.outdir, self.mask_out_file))
        bed_abs_path = os.path.abspath(os.path.join(CWD, self.trf_outdir, f"{self.trf_out_file}.bed"))
        return f"seqtk seq -M {bed_abs_path} {mask_fasta_abs_path} > {self.species}.ssrMasked.fasta"


def argparse_init():
    """初始化命令行参数解析器"""
    parser = argparse.ArgumentParser(description="基因组串联重复序列注释流程")
    parser.add_argument("-g", "--genome", required=True,
                        help="输入基因组FASTA文件路径")
    parser.add_argument("-s", "--species", required=True,
                        help="物种名称（用于输出文件命名）")
    return parser.parse_args()


def main():
    args = argparse_init()
    input_genome = args.genome
    species = args.species
    
    # 检查输入基因组文件
    if not os.path.exists(input_genome):
        print(f"错误：输入基因组文件不存在: {input_genome}")
        sys.exit(1)
    
    # 创建输出目录
    for dirname in [GMATA_OUTDIR, TRF_OUTDIR]:
        folder = File(dirname)
        run_cmd_wait(folder.mkdir())
    
    # 初始化注释流程
    tr_annot = TR_Annotation(input_genome, species)
    
    # 按顺序执行步骤
    run_cmd_wait(tr_annot.run_gmata())
    run_cmd_wait(tr_annot.gmata_mask())
    run_cmd_wait(tr_annot.run_trf())
    run_cmd_wait(tr_annot.run_dta2bed())
    run_cmd_wait(tr_annot.mask_2_lower())
    
    print("\n✅ 所有步骤执行完成！")
    print(f"最终屏蔽基因组文件: {os.path.join(TRF_OUTDIR, f'{species}.ssrMasked.fasta')}")


if __name__ == "__main__":
    # 配置日志
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.FileHandler("tr_annotation.log"), logging.StreamHandler()]
    )
    
    try:
        main()
    except KeyboardInterrupt:
        logging.info("用户中断程序")
        sys.stderr.write("\n用户中断程序！再见！\n")
        sys.exit(0)
    except Exception as e:
        logging.exception("程序执行出错")
        sys.stderr.write(f"\n程序执行出错: {str(e)}\n")
        sys.exit(1)