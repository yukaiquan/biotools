#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed 7 18 15:07:59 2023
This script is used to annotate the genome structure of the species of interest.
多线程并行版本：支持按染色体并行处理
"""
import os
import sys
import argparse
import logging
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import List, Tuple, Optional


# ===================== 软件路径配置 =====================
GMATA = "/publicssd/share/h13713/soft/GMATA/gmata.pl"
GMATA_PATH = "/publicssd/share/h13713/soft/GMATA"
GMATA_CONFIG = "/publicssd/share/h13713/soft/GMATA/default_cfg.txt"

# 自定义GMATA配置（禁用电子PCR以避免内存映射错误）
CUSTOM_GMATA_CONFIG = None  # 将在运行时创建
TRF = "/public/software/apps/geta/biosoft/RepeatMasker-4.1.6/trf"
GMATA_MASK = "/publicssd/share/h13713/soft/GMATA/gssrmsk.pl"
DTA2BED = "/publicssd/share/h13713/soft/dta2bed.py"
SAMTOOLS = "/public/software/apps/geta/biosoft/samtools-1.17/bin/samtools"
BGZIP = "/public/software/apps/geta/biosoft/htslib-1.17/bin/bgzip"
SEQTK = "/publicssd/share/h13713/soft/miniconda3/envs/bwa/bin/seqtk"
RSCRIPT = "/public/software/apps/R/4.4.1/bin/Rscript"

# ===================== 输出目录配置 =====================
GMATA_OUTDIR = "01_GMATA"
TRF_OUTDIR = "02_TRF"

# 自定义GMATA配置文件路径（全局常量）
GMATA_CUSTOM_CONFIG = "gmata_custom.cfg"

# ===================== 全局环境变量设置 =====================
def check_and_update_bashrc():
    """检查并更新 .bashrc 文件中的 GMATA PATH 设置"""
    bashrc_path = os.path.expanduser("~/.bashrc")
    gmata_path_export = f'export PATH="/publicssd/share/h13713/soft/GMATA:$PATH"'
    
    # 如果 .bashrc 文件不存在，创建它
    if not os.path.exists(bashrc_path):
        logging.info(f"创建 .bashrc 文件: {bashrc_path}")
        with open(bashrc_path, 'w') as f:
            f.write(f"# GMATA 环境变量\n{gmata_path_export}\n")
        print(f"\n已在 {bashrc_path} 中添加 GMATA 路径")
        print("请运行 'source ~/.bashrc' 或重新登录以使更改生效")
        return
    
    # 读取 .bashrc 内容
    with open(bashrc_path, 'r') as f:
        bashrc_content = f.read()
    
    # 检查是否已包含 GMATA PATH 设置
    if gmata_path_export in bashrc_content:
        logging.info(".bashrc 中已包含 GMATA PATH 设置")
        return
    
    # 如果没有，追加到文件末尾
    logging.info(f"在 .bashrc 中添加 GMATA PATH 设置")
    with open(bashrc_path, 'a') as f:
        f.write(f"\n# GMATA 环境变量 (自动添加于 {os.popen('date +%Y-%m-%d').read().strip()})\n")
        f.write(f"{gmata_path_export}\n")
    
    print(f"\n已在 {bashrc_path} 中添加 GMATA 路径")
    print("请运行 'source ~/.bashrc' 或重新登录以使更改生效")

# 运行时检查并更新 .bashrc
check_and_update_bashrc()

# 将所有工具路径添加到系统PATH
tool_paths = [
    GMATA_PATH,  # GMATA 路径
    os.path.dirname(SAMTOOLS),  # samtools 路径
    os.path.dirname(BGZIP),  # bgzip 路径
    os.path.dirname(SEQTK),  # seqtk 路径
    os.path.dirname(RSCRIPT),  # Rscript 路径
    os.path.dirname(TRF),  # TRF 路径
]

# 合并所有路径到 PATH 环境变量
for tool_path in tool_paths:
    os.environ["PATH"] = f"{tool_path}:{os.environ.get('PATH', '')}"

# 将GMATA路径添加到Perl模块搜索路径（解决"Can't locate XXX.pm"错误）
os.environ["PERL5LIB"] = f"{GMATA_PATH}:{os.environ.get('PERL5LIB', '')}"

CWD = os.getcwd()


def create_gmata_config() -> str:
    """创建自定义GMATA配置文件，禁用电子PCR以避免内存映射错误"""
    config_content = """[set]:gmat
# GMATA配置文件 - 自动生成
# 禁用电子PCR模块以避免并行运行时的内存映射错误
ModulRun = Y
-r = 5
-m = 2
-x = 10
-s = 0

[set]:ssrfig
# 保留SSR图形化模块（可能影响.fms文件生成）
ModulRun = N

[set]:ssr2gff
ModulRun = N

[set]:ssr2gff3
ModulRun = N

[set]:doprimer_smt
# 保留引物设计模块（可能影响后续步骤）
ModulRun = Y
-p = w
-fl = 200
-l = 2000
-smin = 120
-Smax = 400
-tm = 60
-dir=/public/share/acfaa2ssz7/miniconda3/envs/RepeatMasker/bin/primer3_core
-ts = F
-n = MK

[set]:elctPCR
# 禁用电子PCR定位（导致内存映射错误）
ModulRun = N
-i = 
-o = eMap
-pf = w
-w = 12
-f = 1
-m = 3000
-d = 100-1000
-n = 0
-g = 0
-p = -

[set]:mk2gff3
ModulRun = N
-o = eMap
"""
    
    config_file = os.path.join(CWD, GMATA_CUSTOM_CONFIG)
    with open(config_file, 'w') as f:
        f.write(config_content)
    
    logging.info(f"创建自定义GMATA配置文件: {config_file}")
    return config_file


def get_chromosomes(genome: str) -> List[str]:
    """从基因组索引文件 (.fai) 或使用 samtools 获取染色体名称"""
    fai_file = f"{genome}.fai"
    
    # 如果索引文件存在，直接从索引文件读取（最快）
    if os.path.exists(fai_file):
        chroms = []
        with open(fai_file, 'r') as f:
            for line in f:
                chrom = line.split('\t')[0]  # fai文件格式：chrom\tlength\toffset\tlinebases\tlinewidth
                chroms.append(chrom)
        return chroms
    
    # 如果没有索引文件，使用 samtools faidx 获取染色体列表
    cmd = f"{SAMTOOLS} faidx {genome}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    # 从 .fai 输出中解析染色体名称
    chroms = []
    for line in result.stdout.strip().split('\n'):
        if line:
            chrom = line.split('\t')[0]
            chroms.append(chrom)
    
    return chroms

def check_bgzip_format(gz_file: str) -> bool:
    """检查 .gz 文件是否为 bgzip 格式"""
    # bgzip 格式的文件末尾有特殊的 EOF 标记 (28 bytes)
    # 检查文件最后是否有 BGZF EOF 标记
    try:
        with open(gz_file, 'rb') as f:
            f.seek(-28, 2)  # 从文件末尾向前移动28字节
            eof_marker = f.read(28)
            # BGZF EOF 标记: 1f 8b 08 04 00 00 00 00 00 ff 06 00 42 43 02 00 1b 00 03 00 00 00 00 00 00 00 00 00
            return eof_marker[:18] == b'\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00\x1b\x00'
    except Exception:
        return False


def process_gzip_genome(genome: str) -> str:
    """处理 gzip 压缩的基因组文件，确保是 bgzip 格式"""
    if not genome.endswith('.gz'):
        return genome
    
    logging.info(f"检测到 gzip 压缩文件: {genome}")
    
    # 检查是否为 bgzip 格式
    if check_bgzip_format(genome):
        logging.info("文件已是 bgzip 格式，无需转换")
        return genome
    
    # 需要转换为 bgzip 格式
    logging.info("文件不是 bgzip 格式，开始转换...")
    base_name = genome[:-3]  # 去掉 .gz 后缀
    bgzip_file = f"{base_name}.bgzip.gz"
    
    if os.path.exists(bgzip_file):
        logging.info(f"bgzip 文件已存在，跳过转换: {bgzip_file}")
        return bgzip_file
    
    # 解压并重新用 bgzip 压缩
    # 先解压
    decompress_cmd = f"gunzip -c {genome} > {base_name}.tmp"
    run_cmd_wait(decompress_cmd)
    
    # 用 bgzip 压缩
    compress_cmd = f"{BGZIP} -c {base_name}.tmp > {bgzip_file}"
    run_cmd_wait(compress_cmd)
    
    # 删除临时文件
    os.remove(f"{base_name}.tmp")
    
    logging.info(f"转换完成: {bgzip_file}")
    return bgzip_file

def create_fai_index(genome: str) -> None:
    """创建基因组索引文件 (.fai)"""
    fai_file = f"{genome}.fai"
    if os.path.exists(fai_file):
        logging.info(f"索引文件已存在，跳过创建: {fai_file}")
        return
    
    logging.info(f"创建基因组索引文件: {fai_file}")
    cmd = f"{SAMTOOLS} faidx {genome}"
    run_cmd_wait(cmd)


def extract_chromosome(genome: str, chrom: str, output_dir: str) -> str:
    """使用samtools提取单个染色体序列"""
    output_file = os.path.join(output_dir, f"{chrom}.fa")
    if os.path.exists(output_file):
        logging.info(f"染色体序列已存在，跳过提取: {chrom}")
        return output_file
    
    cmd = f"{SAMTOOLS} faidx {genome} {chrom} > {output_file}"
    run_cmd_wait(cmd)
    return output_file


def convert_trf_dat_to_bed(dat_file: str, bed_file: str, chrom: Optional[str] = None) -> None:
    """
    将TRF的dat格式转换为BED格式（0-based坐标）
    
    TRF使用1-based坐标（第1个碱基是位置1）
    BED使用0-based half-open坐标（第1个碱基是位置0）
    
    转换规则：
    - bed_start = trf_start - 1
    - bed_end = trf_end
    
    TRF dat格式示例：
    27589 27646 18 2 34 10 57 87.0 15 0 0 18 - - - - -AGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAG
    
    字段说明：
    1. Start: 27589 (1-based)
    2. End: 27646 (1-based)
    3-15: 其他TRF参数
    
    BED格式：chrom start end name period size copies consensus
    """
    logging.info(f"转换TRF dat到BED格式: {dat_file} -> {bed_file}")
    
    if not os.path.exists(dat_file):
        logging.error(f"TRF dat文件不存在: {dat_file}")
        return
    
    with open(dat_file, 'r') as f_in, open(bed_file, 'w') as f_out:
        line_num = 0
        for line in f_in:
            line = line.strip()
            
            # 跳过空行和序列行
            if not line or line.startswith('Sequence:'):
                continue
            
            # 跳过参数行（包含多个空格分隔的数字，但没有足够字段的行）
            parts = line.split()
            if len(parts) < 15:
                continue
            
            try:
                # TRF dat格式：Start End Period Size Copies ...
                trf_start = int(parts[0])  # 1-based
                trf_end = int(parts[1])    # 1-based
                period = parts[2]
                size = parts[3]
                copies = parts[4]
                
                # 坐标转换：1-based -> 0-based
                bed_start = trf_start - 1
                bed_end = trf_end
                
                # 使用染色体名称或从文件推断
                chrom_name = chrom if chrom else "chr1"
                
                # BED格式：chrom start end name (可选字段)
                # name格式：period={period};size={size};copies={copies}
                name = f"period={period};size={size};copies={copies}"
                
                f_out.write(f"{chrom_name}\t{bed_start}\t{bed_end}\t{name}\n")
                line_num += 1
            except (ValueError, IndexError) as e:
                logging.warning(f"跳过无效行: {line} (错误: {e})")
                continue
    
    logging.info(f"转换完成，共 {line_num} 条记录")


def mask_ssr_regions(fasta_file: str, ssr_file: str, output_file: str) -> None:
    """
    根据SSR结果文件屏蔽基因组中的SSR区域（转换为小写）
    
    Args:
        fasta_file: 输入FASTA文件
        ssr_file: GMATA的SSR结果文件（.ssr格式）
        output_file: 输出FASTA文件
    """
    logging.info(f"屏蔽SSR区域: {fasta_file} -> {output_file}")
    
    if not os.path.exists(ssr_file):
        logging.error(f"SSR文件不存在: {ssr_file}")
        return
    
    # 解析SSR文件获取位置信息
    ssr_regions = []
    with open(ssr_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                try:
                    # SSR文件格式：chrom start end ...
                    chrom = parts[0]
                    start = int(parts[1])  # 1-based
                    end = int(parts[2])    # 1-based
                    ssr_regions.append((chrom, start - 1, end))  # 转换为0-based
                except (ValueError, IndexError):
                    continue
    
    logging.info(f"读取到 {len(ssr_regions)} 个SSR区域")
    
    # 读取FASTA序列并屏蔽SSR区域
    from collections import defaultdict
    sequences = defaultdict(list)
    current_chrom = None
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>'):
                current_chrom = line[1:].split()[0]
                sequences[current_chrom].append(line)
            else:
                sequences[current_chrom].append(line)
    
    # 将序列合并为字符串
    for chrom in sequences:
        header = sequences[chrom][0]
        seq = ''.join(sequences[chrom][1:])
        sequences[chrom] = (header, seq)
    
    # 屏蔽SSR区域
    for chrom, start, end in ssr_regions:
        if chrom in sequences:
            header, seq = sequences[chrom]
            # 将SSR区域转换为小写
            seq = seq[:start] + seq[start:end].lower() + seq[end:]
            sequences[chrom] = (header, seq)
    
    # 写入输出文件
    with open(output_file, 'w') as f:
        for chrom in sequences:
            header, seq = sequences[chrom]
            f.write(header + '\n')
            # 每行60个字符
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + '\n')
    
    logging.info(f"SSR屏蔽完成: {output_file}")


def run_cmd_wait(cmd: str) -> int:
    """运行命令并等待完成，检查返回码"""
    if cmd == "NA":
        logging.info("跳过已完成的步骤")
        return 0
    
    logging.info(f"执行命令: {cmd}")
    print(f"\n[执行] {cmd}")
    
    # TRF特殊处理：TRF即使成功也可能返回退出码1
    is_trf_cmd = "trf" in cmd.lower() and "2 7 7 80 10 50 500" in cmd
    
    # 构建完整的环境变量（确保子进程能找到所有工具）
    env = os.environ.copy()
    
    # 添加 GMATA 路径到 PATH 和 PERL5LIB（解决并行子进程环境变量问题）
    env["PATH"] = f"{GMATA_PATH}:{env.get('PATH', '')}"
    env["PERL5LIB"] = f"{GMATA_PATH}:{env.get('PERL5LIB', '')}"
    
    try:
        result = subprocess.run(
            cmd, 
            shell=True, 
            check=not is_trf_cmd,  # TRF命令不使用check=True
            capture_output=True,
            text=True,
            env=env  # 传递完整的环境变量
        )
        
        if result.stdout:
            logging.info(f"标准输出:\n{result.stdout}")
        if result.stderr:
            logging.info(f"标准错误:\n{result.stderr}")
        
        # 对于TRF命令，检查输出文件是否存在来判断是否成功
        if is_trf_cmd:
            # TRF的输出信息在stderr中，显示"Done"表示成功
            if "Done" in result.stderr or result.returncode == 0:
                logging.info("TRF执行成功")
                return 0
            else:
                logging.error(f"TRF执行失败！返回码: {result.returncode}")
                logging.error(f"错误输出:\n{result.stderr}")
                print(f"\n[错误] TRF执行失败: {cmd}")
                print(f"错误信息: {result.stderr}")
                sys.exit(1)
        
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
    """串联重复序列注释类 - 支持单染色体处理"""
    def __init__(self, genome: str, species: str, chrom: Optional[str] = None):
        # 每次初始化回到运行目录
        os.chdir(CWD)
        self.genome = genome
        self.genome_basename = os.path.basename(genome)
        self.chrom = chrom
        self.outdir = GMATA_OUTDIR
        self.trf_outdir = TRF_OUTDIR
        self.species = species
        
        # 如果是单染色体模式，使用染色体命名
        if chrom:
            # 先设置染色体fasta文件路径
            self.chrom_fasta = os.path.join(CWD, self.outdir, f"{chrom}.fa")
            # GMATA 输出文件名基于输入文件名：chr1A.fa → chr1A.fa.ssr
            chrom_fa_basename = os.path.basename(self.chrom_fasta)  # chr1A.fa
            self.gmata_out_file = f"{chrom_fa_basename}.ssr"  # chr1A.fa.ssr
            self.mask_out_file = f"{chrom_fa_basename}.ssrMasked.fasta"  # chr1A.fa.ssrMasked.fasta
            self.trf_out_file = f"{self.mask_out_file}.2.7.7.80.10.50.500.dat"
        else:
            self.gmata_out_file = f"{self.genome_basename}.ssr"
            self.mask_out_file = f"{self.genome_basename}.ssrMasked.fasta"
            self.trf_out_file = f"{self.mask_out_file}.2.7.7.80.10.50.500.dat"
            self.chrom_fasta = os.path.join(CWD, genome)

    def run_gmata(self) -> str:
        """运行GMATA进行SSR预测"""
        os.chdir(os.path.join(CWD, self.outdir))
        output_path = os.path.join(CWD, self.outdir, self.gmata_out_file)
        
        if os.path.exists(output_path):
            logging.info(f"GMATA结果已存在，跳过: {output_path}")
            return "NA"
        
        # 使用染色体fasta文件（单染色体模式）或原始基因组文件
        genome_abs_path = os.path.abspath(self.chrom_fasta)
        
        # 使用默认配置文件
        config_file = GMATA_CONFIG
        
        # 将 GMATA 目录添加到 PATH，使其内部的 system() 调用能找到其他脚本
        return f"export PATH={GMATA_PATH}:$PATH && perl {GMATA} -c {config_file} -i {genome_abs_path}"

    def gmata_mask(self) -> str:
        """使用GMATA结果屏蔽基因组中的SSR序列"""
        os.chdir(os.path.join(CWD, self.outdir))
        output_path = os.path.join(CWD, self.outdir, self.mask_out_file)
        
        if os.path.exists(output_path):
            logging.info(f"SSR屏蔽结果已存在，跳过: {output_path}")
            return "NA"
        
        # 首先尝试使用GMATA的gssrmsk.pl
        genome_abs_path = os.path.abspath(self.chrom_fasta)
        ssr_abs_path = os.path.abspath(os.path.join(CWD, self.outdir, self.gmata_out_file))
        
        # 检查.fms文件是否存在（gssrmsk.pl需要）
        fms_file = genome_abs_path + ".fms"
        if os.path.exists(fms_file):
            # 如果.fms文件存在，使用GMATA的屏蔽工具
            return f"export PATH={GMATA_PATH}:$PATH && perl {GMATA_MASK} -i {genome_abs_path} -s 1 -r {ssr_abs_path}"
        else:
            # 否则使用Python实现的SSR屏蔽
            logging.info("未找到.fms文件，使用Python实现的SSR屏蔽")
            output_abs_path = os.path.abspath(output_path)
            mask_ssr_regions(genome_abs_path, ssr_abs_path, output_abs_path)
            return "NA"

    def run_trf(self) -> str:
        """运行TRF进行串联重复序列预测"""
        os.chdir(os.path.join(CWD, self.trf_outdir))
        output_path = os.path.join(CWD, self.trf_outdir, self.trf_out_file)
        
        if os.path.exists(output_path):
            logging.info(f"TRF结果已存在，跳过: {output_path}")
            return "NA"
        
        mask_fasta_abs_path = os.path.abspath(os.path.join(CWD, self.outdir, self.mask_out_file))
        return f"{TRF} {mask_fasta_abs_path} 2 7 7 80 10 50 500 -d -h"

    def process_single_chrom(self) -> Tuple[Optional[str], bool]:
        """处理单个染色体的完整流程（用于并行）"""
        if self.chrom is None:
            logging.error("染色体名称不能为空")
            return (None, False)
        
        try:
            logging.info(f"开始处理染色体: {self.chrom}")
            run_cmd_wait(self.run_gmata())
            run_cmd_wait(self.gmata_mask())
            run_cmd_wait(self.run_trf())
            run_cmd_wait(self.run_dta2bed())
            run_cmd_wait(self.mask_2_lower())
            logging.info(f"染色体 {self.chrom} 处理完成")
            return (self.chrom, True)
        except Exception as e:
            logging.error(f"染色体 {self.chrom} 处理失败: {str(e)}")
            return (self.chrom, False)

    def run_dta2bed(self) -> str:
        """将TRF的dat格式转换为bed格式（0-based坐标）"""
        os.chdir(os.path.join(CWD, self.trf_outdir))
        output_path = os.path.join(CWD, self.trf_outdir, f"{self.trf_out_file}.bed")
        
        if os.path.exists(output_path):
            logging.info(f"BED文件已存在，跳过: {output_path}")
            return "NA"
        
        dat_path = os.path.join(CWD, self.trf_outdir, self.trf_out_file)
        
        # 直接调用转换函数
        convert_trf_dat_to_bed(dat_path, output_path, self.chrom)
        return "NA"  # 已经处理完成，返回NA跳过命令执行

    def mask_2_lower(self) -> str:
        """使用seqtk将bed区域转换为小写字母"""
        os.chdir(os.path.join(CWD, self.trf_outdir))
        # 最终输出文件名：去掉 .fa 后缀，使用染色体名称
        if self.chrom:
            output_name = f"{self.chrom}.ssrMasked.fasta"
        else:
            output_name = f"{self.species}.ssrMasked.fasta"
        
        output_path = os.path.join(CWD, self.trf_outdir, output_name)
        
        if os.path.exists(output_path):
            logging.info(f"最终屏蔽基因组已存在，跳过: {output_path}")
            return "NA"
        
        mask_fasta_abs_path = os.path.abspath(os.path.join(CWD, self.outdir, self.mask_out_file))
        bed_abs_path = os.path.abspath(os.path.join(CWD, self.trf_outdir, f"{self.trf_out_file}.bed"))
        return f"{SEQTK} seq -M {bed_abs_path} {mask_fasta_abs_path} > {output_name}"


def argparse_init():
    """初始化命令行参数解析器"""
    parser = argparse.ArgumentParser(description="基因组串联重复序列注释流程（支持多线程）")
    parser.add_argument("-g", "--genome", required=True,
                        help="输入基因组FASTA文件路径")
    parser.add_argument("-s", "--species", required=True,
                        help="物种名称（用于输出文件命名）")
    parser.add_argument("-t", "--threads", type=int, default=4,
                        help="并行线程数（默认: 4）")
    parser.add_argument("--no-parallel", action="store_true",
                        help="禁用并行处理（使用单线程模式）")
    return parser.parse_args()


def merge_chrom_results(chroms: List[str], species: str) -> None:
    """合并所有染色体的结果文件"""
    logging.info("开始合并各染色体结果...")
    
    # 合并SSR结果
    ssr_out = os.path.join(GMATA_OUTDIR, f"{species}.ssr")
    with open(ssr_out, 'w') as out_f:
        header_written = False
        for chrom in chroms:
            # 文件名：chr1A.fa.ssr
            ssr_file = os.path.join(GMATA_OUTDIR, f"{chrom}.fa.ssr")
            if not os.path.exists(ssr_file):
                logging.warning(f"SSR文件不存在: {ssr_file}")
                continue
            
            with open(ssr_file, 'r') as in_f:
                for line in in_f:
                    if line.startswith('#') and not header_written:
                        out_f.write(line)
                        header_written = True
                    elif not line.startswith('#'):
                        out_f.write(line)
    
    # 合并屏蔽后的基因组文件
    masked_fasta = os.path.join(TRF_OUTDIR, f"{species}.ssrMasked.fasta")
    with open(masked_fasta, 'w') as out_f:
        for chrom in chroms:
            # 文件名：chr1A.ssrMasked.fasta
            chrom_fasta = os.path.join(TRF_OUTDIR, f"{chrom}.ssrMasked.fasta")
            if not os.path.exists(chrom_fasta):
                logging.warning(f"屏蔽基因组文件不存在: {chrom_fasta}")
                continue
            
            with open(chrom_fasta, 'r') as in_f:
                out_f.write(in_f.read())
    
    # 合并BED文件
    bed_out = os.path.join(TRF_OUTDIR, f"{species}.trf.bed")
    with open(bed_out, 'w') as out_f:
        for chrom in chroms:
            # 文件名：chr1A.fa.ssrMasked.fasta.2.7.7.80.10.50.500.dat.bed
            bed_file = os.path.join(TRF_OUTDIR, f"{chrom}.fa.ssrMasked.fasta.2.7.7.80.10.50.500.dat.bed")
            if not os.path.exists(bed_file):
                logging.warning(f"BED文件不存在: {bed_file}")
                continue
            
            with open(bed_file, 'r') as in_f:
                out_f.write(in_f.read())
    
    logging.info("结果合并完成！")


def main():
    args = argparse_init()
    input_genome = args.genome
    species = args.species
    threads = args.threads
    use_parallel = not args.no_parallel
    
    # 使用默认GMATA配置
    
    # 检查输入基因组文件
    if not os.path.exists(input_genome):
        print(f"错误：输入基因组文件不存在: {input_genome}")
        sys.exit(1)
    
    # 处理 gzip 压缩文件（如果是 gzip 格式，确保是 bgzip 格式）
    processed_genome = process_gzip_genome(input_genome)
    
    # 创建输出目录
    for dirname in [GMATA_OUTDIR, TRF_OUTDIR]:
        folder = File(dirname)
        run_cmd_wait(folder.mkdir())
    
    # 创建基因组索引文件（并行处理前必需）
    create_fai_index(processed_genome)
    
    # 提取所有染色体
    chroms = get_chromosomes(processed_genome)
    logging.info(f"检测到 {len(chroms)} 条染色体")
    
    if use_parallel and len(chroms) > 1:
        # ===== 多线程并行模式 =====
        logging.info(f"使用 {threads} 个线程并行处理...")
        
        # 并行提取所有染色体序列
        logging.info(f"并行提取染色体序列（使用 {threads} 个线程）...")
        with ProcessPoolExecutor(max_workers=threads) as executor:
            extract_futures = []
            for chrom in chroms:
                future = executor.submit(extract_chromosome, processed_genome, chrom, GMATA_OUTDIR)
                extract_futures.append((future, chrom))
            
            # 等待所有提取完成并显示进度
            extracted_count = 0
            for future, chrom in extract_futures:
                try:
                    future.result()
                    extracted_count += 1
                    print(f"  [{extracted_count}/{len(chroms)}] 提取完成: {chrom}")
                except Exception as e:
                    logging.error(f"提取染色体 {chrom} 失败: {str(e)}")
        
        logging.info(f"染色体序列提取完成，共 {extracted_count} 条")
        
        # 并行处理每个染色体（GMATA + TRF 分析）
        logging.info("并行执行串联重复序列分析...")
        results = []
        with ProcessPoolExecutor(max_workers=threads) as executor:
            futures = {}
            for chrom in chroms:
                tr_annot = TR_Annotation(processed_genome, species, chrom)
                future = executor.submit(tr_annot.process_single_chrom)
                futures[future] = chrom
            
            for future in as_completed(futures):
                chrom, success = future.result()
                if chrom is not None:  # 只记录有效的染色体结果
                    results.append((chrom, success))
                    if success:
                        print(f"[成功] {chrom} 完成")
                    else:
                        print(f"[失败] {chrom} 失败")
        
        # 合并结果
        successful_chroms = [chrom for chrom, success in results if success and chrom is not None]
        if successful_chroms:
            merge_chrom_results(successful_chroms, species)
        
        # 统计结果
        failed_count = len(chroms) - len(successful_chroms)
        print(f"\n处理统计：")
        print(f"  成功: {len(successful_chroms)} 条染色体")
        print(f"  失败: {failed_count} 条染色体")
    
    else:
        # ===== 单线程模式 =====
        logging.info("使用单线程模式处理...")
        tr_annot = TR_Annotation(processed_genome, species)
        
        # 按顺序执行步骤
        run_cmd_wait(tr_annot.run_gmata())
        run_cmd_wait(tr_annot.gmata_mask())
        run_cmd_wait(tr_annot.run_trf())
        run_cmd_wait(tr_annot.run_dta2bed())
        run_cmd_wait(tr_annot.mask_2_lower())
    
    print("\n所有步骤执行完成！")
    if use_parallel and len(chroms) > 1:
        print(f"合并SSR结果: {os.path.join(GMATA_OUTDIR, f'{species}.ssr')}")
        print(f"合并屏蔽基因组: {os.path.join(TRF_OUTDIR, f'{species}.ssrMasked.fasta')}")
        print(f"合并BED文件: {os.path.join(TRF_OUTDIR, f'{species}.trf.bed')}")
    else:
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