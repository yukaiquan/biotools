#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NGS 测序数据比对流程 (支持单端/双端)
功能: Fastp质控 -> Minimap2比对 -> Samtools排序 -> [Sambamba去重] -> BAM转CRAM
"""
import sys
import os
import shutil
from subprocess import CalledProcessError
import argparse
import math

# ====================== 工具路径配置 ======================
SAMTOOLS = '/publicssd/share/h13713/soft/miniconda3/envs/bwa/bin/samtools'
FASTP = '/publicssd/share/h13713/soft/miniconda3/envs/bwa/bin/fastp'
MINIMAP2 = '/publicssd/share/h13713/soft/miniconda3/envs/bwa/bin/minimap2'
SAMBAMBA = '/publicssd/share/h13713/soft/miniconda3/envs/bwa/bin/sambamba'

# ====================== 工具函数 ======================
def check_file(file_path: str, check_size: bool = True) -> None:
    if not os.path.exists(file_path):
        print(f'Error: File not exists - {file_path}', file=sys.stderr)
        sys.exit(1)
    if check_size and os.path.getsize(file_path) == 0:
        print(f'Error: Empty file - {file_path}', file=sys.stderr)
        sys.exit(1)

def check_dir(dir_path: str) -> None:
    os.makedirs(dir_path, exist_ok=True)

def rm_file(file_path: str) -> None:
    if os.path.exists(file_path):
        try:
            os.remove(file_path)
        except Exception as e:
            print(f"Warning: Failed to remove {file_path}", file=sys.stderr)

def rm_dir(dir_path: str) -> None:
    if os.path.exists(dir_path):
        try:
            shutil.rmtree(dir_path)
        except Exception as e:
            print(f"Warning: Failed to remove {dir_path}", file=sys.stderr)

def is_file_complete(file_path: str) -> bool:
    return os.path.exists(file_path) and os.path.getsize(file_path) > 0

# ====================== Fastp类 ======================
class Fastp:
    def __init__(self, input1: str, input2: str, output_dir: str, prefix: str, thread: int, is_paired: bool = True):
        self.input1 = input1
        self.input2 = input2
        self.output_dir = output_dir
        self.prefix = prefix
        self.thread = thread
        self.is_paired = is_paired

        self.html_report = f"{output_dir}/{prefix}_fastp.html"
        self.json_report = f"{output_dir}/{prefix}_fastp.json"
        self.clean_1 = f"{output_dir}/{prefix}_clean_1.fq.gz"
        self.clean_2 = f"{output_dir}/{prefix}_clean_2.fq.gz" if is_paired else None
        self.single_clean = f"{output_dir}/{prefix}_clean.fq.gz" if not is_paired else None

    def build_command(self) -> str:
        if self.is_paired:
            if is_file_complete(self.html_report) and is_file_complete(self.clean_1) and is_file_complete(self.clean_2):
                print(f"Skip Fastp (paired): Output files exist")
                return 'NA'
            cmd = f"{FASTP} -i {self.input1} -I {self.input2} -o {self.clean_1} -O {self.clean_2} --html {self.html_report} --json {self.json_report} --thread {self.thread}"
        else:
            if is_file_complete(self.html_report) and is_file_complete(self.single_clean):
                print(f"Skip Fastp (single): Output files exist")
                return 'NA'
            cmd = f"{FASTP} -i {self.input1} -o {self.single_clean} --html {self.html_report} --json {self.json_report} --thread {self.thread}"
        return cmd

# ====================== Minimap2类 ======================
class Minimap2AndSort:
    def __init__(self, reference: str, thread: int, rg_info: str, output_dir: str, prefix: str, is_paired: bool = True):
        self.reference = reference
        self.thread = thread
        self.rg_info = rg_info
        self.output_dir = output_dir
        self.prefix = prefix
        self.is_paired = is_paired

        self.sorted_bam = f"{output_dir}/{prefix}.sorted.bam"
        self.sorted_bam_csi = f"{self.sorted_bam}.csi"
        self.log_file = f"{output_dir}/{prefix}_minimap2.log"

        self.samtools_sort_threads = max(1, math.floor(self.thread * 0.08))
        self.minimap2_threads = self.thread - self.samtools_sort_threads
        self.samtools_sort_mem_per_thread = "50G"

    def build_command(self, input_fastq: str, input_fastq2: str = None) -> str:
        if is_file_complete(self.sorted_bam) and is_file_complete(self.sorted_bam_csi) and is_file_complete(self.log_file):
            return 'NA'

        align_mode = "-ax sr"
        
        if self.is_paired:
            input_files = f"{input_fastq} {input_fastq2}"
        else:
            input_files = input_fastq

        cmd = f"{MINIMAP2} -t {self.minimap2_threads} -I 16G {align_mode} -R \"{self.rg_info}\" {self.reference} {input_files} 2> {self.log_file} | {SAMTOOLS} sort --threads {self.samtools_sort_threads} -m {self.samtools_sort_mem_per_thread} -O BAM -o {self.sorted_bam} - && {SAMTOOLS} index -c -@ {self.thread} {self.sorted_bam}"
        return cmd

# ====================== Sambamba Markdup类 ======================
class SambambaMarkdup:
    def __init__(self, thread: int, tmp_dir: str, output_dir: str, prefix: str):
        self.thread = thread
        self.tmp_dir = tmp_dir
        self.output_dir = output_dir
        self.prefix = prefix
        self.rmdup_bam = f"{output_dir}/{prefix}.rmdup.bam"
        self.rmdup_bam_bai = f"{self.rmdup_bam}.bai"

    def build_command(self, input_bam: str) -> str:
        if is_file_complete(self.rmdup_bam) and is_file_complete(self.rmdup_bam_bai):
            return 'NA'
        check_dir(self.tmp_dir)
        return f"{SAMBAMBA} markdup -t {self.thread} --tmpdir={self.tmp_dir} {input_bam} {self.rmdup_bam}"

# ====================== BAM to CRAM类 ======================
class BamToCram:
    def __init__(self, reference: str, thread: int, output_dir: str, prefix: str):
        self.reference = reference
        self.thread = thread
        self.output_dir = output_dir
        self.prefix = prefix
        self.cram = f"{output_dir}/{prefix}_sort.cram"
        self.cram_csi = f"{self.cram}.crai"

    def build_command(self, input_bam: str) -> str:
        if is_file_complete(self.cram) and is_file_complete(self.cram_csi):
            return 'NA'
        return f"{SAMTOOLS} view --threads {self.thread} -T {self.reference} -C -o {self.cram} {input_bam} && {SAMTOOLS} index -@ {self.thread} {self.cram}"

# ====================== 命令执行类 ======================
class CommandExecutor:
    def __init__(self, command: str, dry_run: bool, step_name: str):
        self.command = command
        self.dry_run = dry_run
        self.step_name = step_name

    def run(self) -> None:
        if self.command == 'NA': return
        print(f"\n[{self.step_name}] {self.command}")
        if self.dry_run:
            print(f"[{self.step_name}] Dry run - command not executed")
            return
        try:
            result = os.system(self.command)
            if result != 0:
                raise CalledProcessError(result, self.command)
            print(f"[{self.step_name}] Success")
        except CalledProcessError as e:
            print(f"Error: {self.step_name} failed", file=sys.stderr)
            sys.exit(1)

# ====================== 主函数 ======================
def main():
    parser = argparse.ArgumentParser(description='NGS Pipeline (单端/双端支持, markdup可选)')
    parser.add_argument('-i', '--input1', type=str, required=True, help='R1 fastq.gz (单端模式时为唯一输入)')
    parser.add_argument('-j', '--input2', type=str, default=None, help='R2 fastq.gz (双端模式必填)')
    parser.add_argument('-o', '--output_dir', type=str, required=True)
    parser.add_argument('-r', '--reference', type=str, required=True)
    parser.add_argument('-t', '--thread', type=int, default=8)
    parser.add_argument('-p', '--prefix', type=str, default='sample')
    parser.add_argument('--rg-id', type=str, required=True)
    parser.add_argument('--rg-pl', type=str, default='ILLUMINA')
    parser.add_argument('--rg-lb', type=str, required=True)
    parser.add_argument('--rg-sm', type=str, required=True)
    parser.add_argument('--dry-run', action='store_true')
    parser.add_argument('--keep-intermediate', action='store_true')
    parser.add_argument('--no-markdup', action='store_true', help='跳过 markdup 步骤')

    args = parser.parse_args()

    is_paired = args.input2 is not None
    mode_name = "Paired-end" if is_paired else "Single-end"

    print('#' * 80)
    print(f'  Mode: {mode_name}')
    print(f'  Markdup: {not args.no_markdup}')
    print('#' * 80)

    check_file(args.input1)
    if is_paired:
        check_file(args.input2)
    check_file(args.reference)
    check_file(f"{args.reference}.fai")
    check_dir(args.output_dir)

    # Step 1: Fastp
    print("\n[Step 1] Fastp QC")
    fastp = Fastp(args.input1, args.input2, args.output_dir, args.prefix, args.thread, is_paired)
    CommandExecutor(fastp.build_command(), args.dry_run, "Fastp").run()

    if is_paired:
        clean_1, clean_2 = fastp.clean_1, fastp.clean_2
        if not args.dry_run:
            check_file(clean_1), check_file(clean_2)
        print(f"Fastp output: {clean_1}, {clean_2}")
    else:
        single_clean = fastp.single_clean
        if not args.dry_run:
            check_file(single_clean)
        print(f"Fastp output: {single_clean}")

    # Step 2: Minimap2
    print("\n[Step 2] Minimap2 + Samtools Sort")
    rg_info = f"@RG\\tID:{args.rg_id}\\tLB:{args.rg_lb}\\tPL:{args.rg_pl}\\tPU:{args.rg_id}\\tSM:{args.rg_sm}"
    minimap2_sort = Minimap2AndSort(args.reference, args.thread, rg_info, args.output_dir, args.prefix, is_paired)

    if is_paired:
        cmd = minimap2_sort.build_command(clean_1, clean_2)
    else:
        cmd = minimap2_sort.build_command(single_clean)

    CommandExecutor(cmd, args.dry_run, "Minimap2").run()

    if not args.dry_run:
        check_file(minimap2_sort.sorted_bam), check_file(minimap2_sort.sorted_bam_csi)
    print(f"BAM: {minimap2_sort.sorted_bam}")

    # Step 3: Markdup (可选)
    do_markdup = not args.no_markdup
    if do_markdup:
        print("\n[Step 3] Sambamba Markdup")
        tmp_dir = os.path.join(args.output_dir, "tmp_markdup")
        markdup = SambambaMarkdup(args.thread, tmp_dir, args.output_dir, args.prefix)
        CommandExecutor(markdup.build_command(minimap2_sort.sorted_bam), args.dry_run, "Markdup").run()

        if not args.dry_run:
            check_file(markdup.rmdup_bam), check_file(markdup.rmdup_bam_bai)
        print(f"Rmdup BAM: {markdup.rmdup_bam}")
        final_bam = markdup.rmdup_bam
    else:
        print("\n[Step 3] Skip Markdup (--no-markdup)")
        final_bam = minimap2_sort.sorted_bam

    # Step 4: CRAM
    print("\n[Step 4] BAM -> CRAM")
    bam_to_cram = BamToCram(args.reference, args.thread, args.output_dir, args.prefix)
    CommandExecutor(bam_to_cram.build_command(final_bam), args.dry_run, "BAM->CRAM").run()

    final_cram = bam_to_cram.cram
    if not args.dry_run:
        check_file(final_cram), check_file(f"{final_cram}.crai")
    print(f"Final CRAM: {final_cram}")

    # Cleanup
    if not args.keep_intermediate and not args.dry_run:
        print("\n[Cleanup]")
        if is_paired:
            rm_file(fastp.clean_1)
            rm_file(fastp.clean_2)
        else:
            rm_file(fastp.single_clean)
        rm_file(minimap2_sort.sorted_bam)
        rm_file(minimap2_sort.sorted_bam_csi)
        if do_markdup:
            rm_file(markdup.rmdup_bam)
            rm_file(markdup.rmdup_bam_bai)
            rm_dir(tmp_dir)
        rm_file(minimap2_sort.log_file)
        print("Cleanup done")

    print("\n" + '#' * 80)
    print(f"Pipeline completed! Mode: {mode_name}")
    print(f"Final: {final_cram}")
    print('#' * 80)

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write('\nInterrupted\n')
        sys.exit(1)
    except CalledProcessError as e:
        sys.stderr.write(f'\nError: {e}\n')
        sys.exit(1)
