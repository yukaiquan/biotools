#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NGS Sequencing Data Alignment Pipeline (Supports Single-end/Paired-end)
Functions: [Fastp QC (Optional)] -> Minimap2 Alignment -> Samtools Sort -> [Sambamba Markdup] -> BAM to CRAM (Optional)
New: --uniq-mapping option to filter non-unique mapping reads
New: --skip-fastp option to skip Fastp QC step
"""
import sys
import os
import shutil
from subprocess import CalledProcessError
import argparse
import math
from typing import Optional

# ====================== Tool Path Configuration ======================
SAMTOOLS = '/publicssd/share/h13713/soft/miniconda3/envs/bwa/bin/samtools'
FASTP = '/publicssd/share/h13713/soft/miniconda3/envs/bwa/bin/fastp'
MINIMAP2 = '/publicssd/share/h13713/soft/miniconda3/envs/bwa/bin/minimap2'
SAMBAMBA = '/publicssd/share/h13713/soft/miniconda3/envs/bwa/bin/sambamba'

# ====================== Utility Functions ======================
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

# ====================== Fastp Class ======================
class Fastp:
    def __init__(self, input1: str, input2: Optional[str], output_dir: str, prefix: str, thread: int, is_paired: bool = True) -> None:
        self.input1 = input1
        self.input2 = input2
        self.output_dir = output_dir
        self.prefix = prefix
        self.thread = thread
        self.is_paired = is_paired

        self.html_report = f"{output_dir}/{prefix}_fastp.html"
        self.json_report = f"{output_dir}/{prefix}_fastp.json"
        self.clean_1: Optional[str] = f"{output_dir}/{prefix}_clean_1.fq.gz"
        self.clean_2: Optional[str] = f"{output_dir}/{prefix}_clean_2.fq.gz" if is_paired else None
        self.single_clean: Optional[str] = f"{output_dir}/{prefix}_clean.fq.gz" if not is_paired else None

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

# ====================== Minimap2 Class ======================
class Minimap2AndSort:
    def __init__(self, reference: str, thread: int, rg_info: str, output_dir: str, prefix: str, 
                 is_paired: bool = True, uniq_mapping: bool = False):
        self.reference = reference
        self.thread = thread
        self.rg_info = rg_info
        self.output_dir = output_dir
        self.prefix = prefix
        self.is_paired = is_paired
        self.uniq_mapping = uniq_mapping

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

        # Basic minimap2 command
        minimap2_cmd = f"{MINIMAP2} -t {self.minimap2_threads} -I 16G {align_mode} -R \"{self.rg_info}\" {self.reference} {input_files} 2> {self.log_file}"
        
        # Unique mapping filter step
        if self.uniq_mapping:
            filter_cmd = f"| {SAMTOOLS} view -h -q 20 -F 3332"
        else:
            filter_cmd = ""
        
        # Sort and index command
        sort_index_cmd = f"| {SAMTOOLS} sort --threads {self.samtools_sort_threads} -m {self.samtools_sort_mem_per_thread} -O BAM -o {self.sorted_bam} - && {SAMTOOLS} index -c -@ {self.thread} {self.sorted_bam}"
        
        # Combine complete command
        cmd = f"{minimap2_cmd} {filter_cmd} {sort_index_cmd}"
        return cmd

# ====================== Sambamba Markdup Class ======================
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

# ====================== BAM to CRAM Class ======================
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

# ====================== Command Executor Class ======================
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

# ====================== Main Function ======================
def main():
    parser = argparse.ArgumentParser(description='NGS Pipeline (Supports single-end/paired-end, markdup optional)')
    parser.add_argument('-i', '--input1', type=str, required=True, help='R1 fastq.gz (single input for single-end mode)')
    parser.add_argument('-j', '--input2', type=str, default=None, help='R2 fastq.gz (required for paired-end mode)')
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
    parser.add_argument('--no-markdup', action='store_true', help='Skip markdup step')
    parser.add_argument('-f', '--output-format', type=str, default='cram', choices=['cram', 'bam'], help='Output format (default: cram, optional: bam)')
    parser.add_argument('--uniq-mapping', action='store_true', help='Enable unique mapping filter (samtools view -q 20 -F 3332)')
    # New: Skip Fastp option
    parser.add_argument('--skip-fastp', action='store_true', help='Skip Fastp QC step and use raw fastq for alignment')

    args = parser.parse_args()

    is_paired = args.input2 is not None
    mode_name = "Paired-end" if is_paired else "Single-end"
    output_format = args.output_format.lower()

    print('#' * 80)
    print(f'  Mode: {mode_name}')
    print(f'  Skip Fastp: {args.skip_fastp}')
    print(f'  Markdup: {not args.no_markdup}')
    print(f'  Output Format: {output_format.upper()}')
    print(f'  Unique Mapping: {args.uniq_mapping}')
    print('#' * 80)

    # Check input files
    check_file(args.input1)
    if is_paired:
        check_file(args.input2)
    check_file(args.reference)
    check_file(f"{args.reference}.fai")
    check_dir(args.output_dir)

    # Define fastq file variables for subsequent steps
    fastq_1: Optional[str] = None
    fastq_2: Optional[str] = None
    fastp_instance: Optional[Fastp] = None

    # Step 1: Fastp (Optional)
    if not args.skip_fastp:
        print("\n[Step 1] Fastp QC")
        fastp_instance = Fastp(args.input1, args.input2, args.output_dir, args.prefix, args.thread, is_paired)
        CommandExecutor(fastp_instance.build_command(), args.dry_run, "Fastp").run()

        if is_paired:
            fastq_1, fastq_2 = fastp_instance.clean_1, fastp_instance.clean_2
            if not args.dry_run:
                assert fastq_1 is not None and fastq_2 is not None, "Fastp output files should not be None"
                check_file(fastq_1), check_file(fastq_2)
            print(f"Fastp output: {fastq_1}, {fastq_2}")
        else:
            fastq_1 = fastp_instance.single_clean
            if not args.dry_run:
                assert fastq_1 is not None, "Fastp output file should not be None"
                check_file(fastq_1)
            print(f"Fastp output: {fastq_1}")
    else:
        print("\n[Step 1] Skip Fastp (--skip-fastp)")
        # Use raw input files directly
        fastq_1 = args.input1
        if is_paired:
            fastq_2 = args.input2
        print(f"Using raw fastq: {fastq_1}" + (f", {fastq_2}" if is_paired else ""))

    # Step 2: Minimap2
    print("\n[Step 2] Minimap2 + Samtools Sort")
    rg_info = f"@RG\\tID:{args.rg_id}\\tLB:{args.rg_lb}\\tPL:{args.rg_pl}\\tPU:{args.rg_id}\\tSM:{args.rg_sm}"
    minimap2_sort = Minimap2AndSort(
        args.reference, args.thread, rg_info, args.output_dir, args.prefix, 
        is_paired, args.uniq_mapping
    )

    assert fastq_1 is not None, "fastq_1 should not be None at Minimap2 step"
    if is_paired:
        assert fastq_2 is not None, "fastq_2 should not be None in paired-end mode"
        cmd = minimap2_sort.build_command(fastq_1, fastq_2)
    else:
        cmd = minimap2_sort.build_command(fastq_1)

    CommandExecutor(cmd, args.dry_run, "Minimap2").run()

    if not args.dry_run:
        check_file(minimap2_sort.sorted_bam), check_file(minimap2_sort.sorted_bam_csi)
    print(f"BAM: {minimap2_sort.sorted_bam}")

    # Step 3: Markdup (Optional)
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

    # Step 4: Convert to CRAM based on output format
    final_file: str = ""
    if output_format == "cram":
        print("\n[Step 4] BAM -> CRAM")
        bam_to_cram = BamToCram(args.reference, args.thread, args.output_dir, args.prefix)
        CommandExecutor(bam_to_cram.build_command(final_bam), args.dry_run, "BAM->CRAM").run()

        final_file = bam_to_cram.cram
        if not args.dry_run:
            check_file(final_file)
            check_file(f"{final_file}.crai")
        print(f"Final CRAM: {final_file}")
    else:
        print("\n[Step 4] Skip BAM->CRAM (Output BAM format)")
        final_file = final_bam
        if not args.dry_run:
            check_file(final_file)
            # Check BAM index
            index_files = [f"{final_file}.bai", f"{final_file}.csi"]
            if not any(os.path.exists(idx) for idx in index_files):
                print(f"[Warning] BAM index not found, creating index...")
                index_cmd = f"{SAMTOOLS} index -@ {args.thread} {final_file}"
                CommandExecutor(index_cmd, args.dry_run, "Samtools Index").run()
        print(f"Final BAM: {final_file}")

    # Cleanup: Only clean intermediate files, not final output files and original input files
    if not args.keep_intermediate and not args.dry_run:
        print("\n[Cleanup]")
        # Only clean clean files if Fastp was executed
        if not args.skip_fastp and fastp_instance is not None:
            if is_paired:
                if fastp_instance.clean_1 is not None:
                    rm_file(fastp_instance.clean_1)
                if fastp_instance.clean_2 is not None:
                    rm_file(fastp_instance.clean_2)
            else:
                if fastp_instance.single_clean is not None:
                    rm_file(fastp_instance.single_clean)
        
        # Clean original sorted BAM (not final file)
        rm_file(minimap2_sort.sorted_bam)
        rm_file(minimap2_sort.sorted_bam_csi)
        
        # If CRAM format, clean intermediate deduplicated BAM; BAM format does not clean final file
        if do_markdup and output_format == "cram":
            rm_file(markdup.rmdup_bam)
            rm_file(markdup.rmdup_bam_bai)
            rm_dir(tmp_dir)  # type: ignore
        
        print("Cleanup done")

    print("\n" + '#' * 80)
    print(f"Pipeline completed! Mode: {mode_name}")
    print(f"Skip Fastp: {args.skip_fastp}")
    print(f"Output Format: {output_format.upper()}")
    print(f"Final File: {final_file}")
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
