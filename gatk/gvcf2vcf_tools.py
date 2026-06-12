#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import os
import shutil
import re
import shlex
from subprocess import CalledProcessError
import argparse
import math
from typing import Optional, List, Tuple, Dict
from pathlib import Path

# ====================== Tool Path Configuration ======================
GATK = '/publicssd/share/h13713/soft/miniconda3/envs/bwa/bin/gatk'
BCFTOOLS = '/publicssd/share/h13713/soft/miniconda3/envs/bwa/bin/bcftools'
BGZIP = '/publicssd/share/h13713/soft/miniconda3/envs/bwa/bin/bgzip'

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
            print(f"Removed: {file_path}")
        except Exception as e:
            print(f"Warning: Failed to remove {file_path}: {e}", file=sys.stderr)

def rm_files(file_list: List[str]) -> None:
    """Remove multiple files"""
    for file_path in file_list:
        rm_file(file_path)

def is_file_complete(file_path: str) -> bool:
    return os.path.exists(file_path) and os.path.getsize(file_path) > 0


def index_vcf(vcf_file: str, dry_run: bool) -> None:
    """
    Create tabix index (.tbi) for a VCF file using bcftools index -t.
    Only creates index if it doesn't exist or is older than the VCF file.
    """
    if not is_file_complete(vcf_file):
        print(f"Warning: skipping index for missing VCF: {vcf_file}")
        return
    
    # Check if index already exists and is up to date
    index_file = f"{vcf_file}.csi"
    if os.path.exists(index_file):
        vcf_mtime = os.path.getmtime(vcf_file)
        index_mtime = os.path.getmtime(index_file)
        if index_mtime >= vcf_mtime:
            print(f"  Index already exists and is up to date: {index_file}")
            return
        else:
            print(f"  Index is outdated, recreating: {index_file}")
            rm_file(index_file)
    
    cmd = f"{BCFTOOLS} index -c {shlex.quote(vcf_file)}"
    executor = CommandExecutor(cmd, dry_run, f"Index {os.path.basename(vcf_file)}")
    executor.run()


def locate_interval_vcf(vcf_dir: str, chrom: str, interval_name: str, vcf_type: str) -> Optional[str]:
    candidates = [
        os.path.join(vcf_dir, chrom, f"{interval_name}.{vcf_type}.vcf.gz"),
        os.path.join(vcf_dir, f"{interval_name}.{vcf_type}.vcf.gz"),
    ]
    for path in candidates:
        if os.path.exists(path):
            return path
    return None


def parse_interval_name_from_vcf(filename: str) -> Optional[str]:
    m = re.match(r'^(?P<interval>.+)\.(?:snp|indel)\.vcf\.gz$', filename)
    if m:
        return m.group('interval')
    return None


def validate_vcf_interval_matches(vcf_dir: str, expected_intervals: set) -> None:
    found_intervals = set()
    extra_intervals = set()
    for root, _, files in os.walk(vcf_dir):
        for filename in files:
            interval_name = parse_interval_name_from_vcf(filename)
            if not interval_name:
                continue
            found_intervals.add(interval_name)
            if interval_name not in expected_intervals:
                extra_intervals.add(interval_name)

    for interval_name in sorted(extra_intervals):
        print(f"Warning: VCF interval not in intervals.txt: {interval_name}")

    missing_intervals = expected_intervals - found_intervals
    for interval_name in sorted(missing_intervals):
        print(f"Warning: expected interval from intervals.txt missing any VCF: {interval_name}")


def infer_chromosome(name: str) -> Optional[str]:
    chromosome_re = re.compile(r'([_\.]?chr[^_\.]+)$', re.IGNORECASE)
    match = chromosome_re.search(name)
    if match:
        return match.group(1).lstrip('_.')
    return None


def parse_gvcf_list(gvcf_list: str) -> List[Tuple[str, str, Optional[str]]]:
    samples = []
    with open(gvcf_list, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                sample_name = parts[0]
                gvcf_path = parts[1]
            elif len(parts) == 1:
                gvcf_path = parts[0]
                sample_name = os.path.basename(gvcf_path).replace('.gvcf.gz', '').replace('.g.vcf.gz', '')
            else:
                continue

            sample_chr = infer_chromosome(sample_name)
            if sample_chr is None:
                filename = os.path.basename(gvcf_path)
                stem = filename
                if stem.endswith('.gvcf.gz'):
                    stem = stem[:-8]
                elif stem.endswith('.g.vcf.gz'):
                    stem = stem[:-9]
                sample_chr = infer_chromosome(stem)

            samples.append((sample_name, gvcf_path, sample_chr))
    return samples


# ====================== Command Executor Class ======================
class CommandExecutor:
    def __init__(self, command: str, dry_run: bool, step_name: str):
        self.command = command
        self.dry_run = dry_run
        self.step_name = step_name

    def run(self) -> None:
        if self.command == 'NA':
            return
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


# ====================== Step 0: Generate Intervals ======================
def generate_intervals(fai_file: str, output_file: str, window_size: int = 5000000, overlap_size: int = 20000) -> None:
    """
    Generate interval file from genome fai file.
    
    Format: chr  region_start  region_end  vcf_start  vcf_end
    - region_start/region_end: used for bcftools view -r (extraction region with overlap)
    - vcf_start/vcf_end: used for final VCF filtering (non-overlapping, to avoid duplicates)
    
    Example with 5Mb window and 20kb overlap:
        chr1H  1        5000000   1        4990000
        chr1H  4980001  10000000  4990001  9990000
        chr1H  9980001  15000000  9990001  14990000
    
    Rationale:
    - Extraction regions overlap (1-5000000, 4980001-10000000) to capture boundary variants
    - VCF filter regions are non-overlapping (1-4990000, 4990001-9990000) to avoid duplicates
    """
    print(f"\n[Step 0] Generating intervals from {fai_file}")
    print(f"  Window size: {window_size:,} bp")
    print(f"  Overlap size: {overlap_size:,} bp")
    
    check_file(fai_file)
    
    intervals = []
    with open(fai_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 2:
                continue
            chrom = parts[0]
            chrom_length = int(parts[1])
            
            # Generate windows
            # Strategy: divide genome at window_size boundaries (5M, 10M, 15M...)
            # - Window N ends at N × window_size
            # - Window N (N>1) starts at (N-1) × window_size - overlap_size + 1
            # - First window starts at 1
            
            window_num = 1
            prev_vcf_end = 0
            
            # Margin for filtering (half of overlap)
            margin = overlap_size // 2  # 10000 for overlap=20000
            
            while True:
                # Calculate window boundaries
                # Window N ends at N × window_size
                window_end = min(window_num * window_size, chrom_length)
                
                if window_num == 1:
                    # First window: starts at 1
                    window_start = 1
                else:
                    # Subsequent windows: start with overlap
                    window_start = (window_num - 1) * window_size - overlap_size + 1
                
                # Break if window_start exceeds chrom_length
                if window_start > chrom_length:
                    break
                
                # Calculate VCF filter coordinates (non-overlapping)
                if window_num == 1:
                    # First window: only remove margin from end
                    vcf_start = window_start
                    vcf_end = window_end - margin
                else:
                    # Subsequent windows: start after previous vcf_end
                    vcf_start = prev_vcf_end + 1
                    vcf_end = window_end - margin
                
                # Ensure vcf_end doesn't exceed chrom_length
                vcf_end = min(vcf_end, chrom_length)
                
                # Save coordinates
                intervals.append((chrom, window_start, window_end, vcf_start, vcf_end))
                
                # Prepare for next window
                prev_vcf_end = vcf_end
                window_num += 1
                
                # Break if we've reached chrom_length
                if window_end >= chrom_length:
                    break
    
    # Write to file
    with open(output_file, 'w') as f:
        for chrom, reg_start, reg_end, vcf_start, vcf_end in intervals:
            f.write(f"{chrom}\t{reg_start}\t{reg_end}\t{vcf_start}\t{vcf_end}\n")
    
    print(f"Generated {len(intervals)} intervals -> {output_file}")


def read_intervals(interval_file: str) -> List[Tuple[str, int, int, int, int]]:
    """Read interval file and return list of intervals"""
    check_file(interval_file)
    
    intervals = []
    with open(interval_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) != 5:
                continue
            chrom = parts[0]
            reg_start = int(parts[1])
            reg_end = int(parts[2])
            vcf_start = int(parts[3])
            vcf_end = int(parts[4])
            intervals.append((chrom, reg_start, reg_end, vcf_start, vcf_end))
    
    return intervals


# ====================== Step 1a: Extract Region GVCFs ======================
def extract_region_gvcfs(
    gvcf_list: str,
    interval_file: str,
    output_dir: str,
    thread: int,
    dry_run: bool,
    print_commands: bool = False
) -> Dict[int, List[str]]:
    """
    Step 1a: Extract GVCF regions using bcftools.
    
    For each interval:
    - If print_commands: generate a command file in interval directory
    - Otherwise: execute bcftools commands
    
    Returns: dict mapping interval_id to list of extracted GVCF paths
    """
    if not print_commands:
        print(f"\n[Step 1a] Extracting GVCF regions")
    
    check_file(gvcf_list)
    check_dir(output_dir)
    
    # Read GVCF list (format: sample_name<TAB>gvcf_path or just gvcf_path)
    samples = parse_gvcf_list(gvcf_list)
    
    if not print_commands:
        print(f"Found {len(samples)} samples")
    
    # Read intervals
    intervals = read_intervals(interval_file)
    if not print_commands:
        print(f"Found {len(intervals)} intervals")
    
    # Process each interval
    interval_gvcfs = {}
    
    for idx, (chrom, reg_start, reg_end, vcf_start, vcf_end) in enumerate(intervals):
        interval_id = idx + 1
        interval_name = f"{chrom}_{reg_start}_{reg_end}"
        
        if not print_commands:
            print(f"\n[Interval {interval_id}] {chrom}:{reg_start}-{reg_end}")
        
        # Create chromosome subdirectory
        chrom_dir = os.path.join(output_dir, chrom)
        check_dir(chrom_dir)
        interval_dir = os.path.join(chrom_dir, interval_name)
        check_dir(interval_dir)
        
        # Extract region from each sample
        region_str = f"{chrom}:{reg_start}-{reg_end}"
        sample_gvcfs = []
        commands = []
        
        # Use chromosome-aware sample selection when possible
        matched_samples = [
            (name, path)
            for name, path, sample_chr in samples
            if sample_chr is None or sample_chr == chrom
        ]

        if not matched_samples:
            print(f"  Warning: no matching sample files found for chromosome {chrom} in gvcf_list")

        for sample_name, gvcf_path in matched_samples:
            output_gvcf = os.path.join(interval_dir, f"{sample_name}.gvcf")

            # Check if already exists
            if is_file_complete(output_gvcf):
                if not print_commands:
                    print(f"  Skip extraction (exists): {sample_name}")
                sample_gvcfs.append(output_gvcf)
                continue

            # Generate bcftools command
            cmd = (
                f"{BCFTOOLS} view --threads {thread} -r {region_str} -Ov "
                f"{shlex.quote(gvcf_path)} -o {shlex.quote(output_gvcf)}"
            )

            if print_commands:
                commands.append(cmd)
            else:
                # Execute command
                executor = CommandExecutor(cmd, dry_run, f"Extract {sample_name}")
                executor.run()

            sample_gvcfs.append(output_gvcf)
        
        # Write commands to file if print_commands mode
        if print_commands and commands:
            cmd_file = os.path.join(interval_dir, "extract_commands.sh")
            with open(cmd_file, 'w') as f:
                for cmd in commands:
                    f.write(cmd + "\n")
            print(f"Generated: {cmd_file} ({len(commands)} commands)")
        
        interval_gvcfs[interval_id] = sample_gvcfs
    
    return interval_gvcfs


# ====================== Step 1b: Combine Region GVCFs ======================
def combine_region_gvcfs(
    gvcf_list: Optional[str],
    interval_file: Optional[str],
    output_dir: Optional[str],
    reference: str,
    java_options: str,
    dry_run: bool,
    keep_intermediate: bool,
    interval_dir: str = None
) -> Dict[str, str]:
    """
    Step 1b: Combine extracted GVCFs for each interval.
    
    Two modes:
    1. Process all intervals: provide gvcf_list, interval_file, output_dir
    2. Process single interval: provide interval_dir, reference
    
    Returns: dict mapping interval region to combined GVCF path
    """
    print(f"\n[Step 1b] Combining GVCFs by interval")
    
    interval_gvcfs = {}
    
    # Mode 1: Process single interval directory
    if interval_dir:
        # Check interval directory first
        if not os.path.isdir(interval_dir):
            print(f"Error: Interval directory not found: {interval_dir}", file=sys.stderr)
            sys.exit(1)
        
        # Then check reference file
        check_file(reference)
        
        print(f"Processing single interval: {interval_dir}")
        
        # Find all .gvcf files in the directory
        sample_gvcfs = []
        for filename in sorted(os.listdir(interval_dir)):
            if filename.endswith('.gvcf'):
                sample_gvcfs.append(os.path.join(interval_dir, filename))
        
        if not sample_gvcfs:
            print(f"Error: No .gvcf files found in {interval_dir}", file=sys.stderr)
            sys.exit(1)
        
        print(f"Found {len(sample_gvcfs)} sample GVCFs")
        
        # Extract interval name from directory (normalize path to handle trailing slash)
        interval_dir = os.path.normpath(interval_dir)
        interval_name = os.path.basename(interval_dir)
        
        # Output to parent directory
        parent_dir = os.path.dirname(interval_dir)
        combined_gvcf = os.path.join(parent_dir, f"{interval_name}.gvcf")
        
        if is_file_complete(combined_gvcf):
            print(f"Skip combine (exists): {combined_gvcf}")
            return {interval_dir: combined_gvcf}
        
        # Build GATK CombineGVCFs command
        vcf_params = ' '.join([f'-V {shlex.quote(gvcf)}' for gvcf in sample_gvcfs])
        cmd = (
            f"{GATK} --java-options \"{java_options}\" CombineGVCFs -R {shlex.quote(reference)} {vcf_params} "
            f"--create-output-variant-index false -O {shlex.quote(combined_gvcf)}"
        )
        
        executor = CommandExecutor(cmd, dry_run, f"CombineGVCFs {interval_name}")
        executor.run()
        
        interval_gvcfs[interval_dir] = combined_gvcf
        
        # Clean intermediate files
        if not keep_intermediate and not dry_run:
            rm_files(sample_gvcfs)
        
        return interval_gvcfs
    
    # Mode 2: Process all intervals
    if gvcf_list is None or interval_file is None or output_dir is None:
        print("Error: Batch mode requires gvcf_list, interval_file, and output_dir", file=sys.stderr)
        sys.exit(1)
    
    check_file(gvcf_list)
    check_file(reference)
    check_dir(output_dir)
    
    # Read GVCF list to get sample names and chromosome annotations
    samples = parse_gvcf_list(gvcf_list)
    print(f"Found {len(samples)} sample entries")
    
    # Read intervals
    intervals = read_intervals(interval_file)
    print(f"Found {len(intervals)} intervals")
    
    # Process each interval
    for idx, (chrom, reg_start, reg_end, vcf_start, vcf_end) in enumerate(intervals):
        interval_id = idx + 1
        interval_name = f"{chrom}_{reg_start}_{reg_end}"
        print(f"\n[Interval {interval_id}] {chrom}:{reg_start}-{reg_end}")
        
        # Create chromosome subdirectory
        chrom_dir = os.path.join(output_dir, chrom)
        interval_dir = os.path.join(chrom_dir, interval_name)
        
        # Only use samples that match this chromosome when chromosome info is available
        matched_samples = [
            (name, path)
            for name, path, sample_chr in samples
            if sample_chr is None or sample_chr == chrom
        ]
        if not matched_samples:
            print(f"  Warning: no matching samples found for chromosome {chrom}; skipping interval")
            continue

        # Check if all matched sample GVCFs exist
        sample_gvcfs = []
        missing_files = []
        for sample_name, _ in matched_samples:
            sample_gvcf = os.path.join(interval_dir, f"{sample_name}.gvcf")
            if not is_file_complete(sample_gvcf):
                missing_files.append(sample_gvcf)
            sample_gvcfs.append(sample_gvcf)
        
        if missing_files:
            print(f"  Warning: {len(missing_files)} sample GVCFs not found. Skipping this interval.")
            for f in missing_files[:5]:  # Show first 5 missing files
                print(f"    - {f}")
            if len(missing_files) > 5:
                print(f"    ... and {len(missing_files) - 5} more")
            continue
        
        # Combine all samples for this interval
        combined_gvcf = os.path.join(chrom_dir, f"{interval_name}.gvcf")
        
        if is_file_complete(combined_gvcf):
            print(f"  Skip combine (exists): {combined_gvcf}")
            interval_gvcfs[f"{chrom}:{reg_start}-{reg_end}"] = combined_gvcf
            
            # Clean intermediate files
            if not keep_intermediate and not dry_run:
                rm_files(sample_gvcfs)
            continue
        
        # Build GATK CombineGVCFs command
        vcf_params = ' '.join([f'-V {shlex.quote(gvcf)}' for gvcf in sample_gvcfs])
        cmd = (
            f"{GATK} --java-options \"{java_options}\" CombineGVCFs -R {shlex.quote(reference)} {vcf_params} "
            f"--create-output-variant-index false -O {shlex.quote(combined_gvcf)}"
        )
        
        executor = CommandExecutor(cmd, dry_run, f"CombineGVCFs interval {interval_id}")
        executor.run()
        
        interval_gvcfs[f"{chrom}:{reg_start}-{reg_end}"] = combined_gvcf
        
        # Clean intermediate files
        if not keep_intermediate and not dry_run:
            rm_files(sample_gvcfs)
    
    return interval_gvcfs


# ====================== Step 2: Genotype and Generate VCF ======================
def genotype_and_generate_vcf(
    gvcf_file: str,
    output_dir: str,
    reference: str,
    thread: int,
    java_options: str,
    dry_run: bool,
    keep_intermediate: bool,
    output_prefix: str = None
) -> Tuple[str, str]:
    """
    Step 2: Genotype a single GVCF file and generate SNP/INDEL VCFs.

    1. Run GenotypeGVCFs on a single GVCF
    2. Separate SNPs and INDELs

    Returns: (snp_vcf_path, indel_vcf_path)
    """
    print(f"\n[Step 2] Genotyping GVCF and generating VCFs")

    check_file(gvcf_file)
    check_file(reference)
    check_dir(output_dir)

    # Extract base name
    gvcf_basename = os.path.basename(gvcf_file).replace('.gvcf', '')

    # Use basename as default prefix if not specified
    if output_prefix is None:
        output_prefix = gvcf_basename

    # Genotype
    output_vcf = os.path.join(output_dir, f"{gvcf_basename}.vcf")

    if is_file_complete(output_vcf):
        print(f"  Skip genotyping (exists): {output_vcf}")
    else:
        cmd = f"{GATK} --java-options \"{java_options}\" GenotypeGVCFs -R {shlex.quote(reference)} -V {shlex.quote(gvcf_file)} --create-output-variant-index false -O {shlex.quote(output_vcf)}"
        executor = CommandExecutor(cmd, dry_run, "GenotypeGVCFs")
        executor.run()

    # Separate SNPs and INDELs
    snp_vcf = os.path.join(output_dir, f"{output_prefix}.snp.vcf.gz")
    indel_vcf = os.path.join(output_dir, f"{output_prefix}.indel.vcf.gz")

    if not is_file_complete(snp_vcf):
        cmd = f"{GATK} --java-options \"{java_options}\" SelectVariants -R {shlex.quote(reference)} -V {shlex.quote(output_vcf)} -select-type SNP --create-output-variant-index false -O {shlex.quote(snp_vcf)}"
        executor = CommandExecutor(cmd, dry_run, "SelectVariants SNP")
        executor.run()

    if not is_file_complete(indel_vcf):
        cmd = f"{GATK} --java-options \"{java_options}\" SelectVariants -R {shlex.quote(reference)} -V {shlex.quote(output_vcf)} -select-type INDEL --create-output-variant-index false -O {shlex.quote(indel_vcf)}"
        executor = CommandExecutor(cmd, dry_run, "SelectVariants INDEL")
        executor.run()

    # Auto-index the output SNP and INDEL VCFs
    index_vcf(snp_vcf, dry_run)
    index_vcf(indel_vcf, dry_run)

    # Clean intermediate VCF
    if not keep_intermediate and not dry_run:
        rm_file(output_vcf)

    return snp_vcf, indel_vcf


# ====================== Step 3: Merge Interval VCFs ======================
def merge_interval_vcfs(
    interval_file: str,
    vcf_dir: str,
    output_dir: str,
    snp_header: str,
    indel_header: str,
    thread: int,
    dry_run: bool
) -> Dict[str, Tuple[str, str]]:
    """
    Step 3: Merge interval VCFs by chromosome using vcf_start/vcf_end filtering.
    Supports separate header files for SNP and INDEL.
    
    For each chromosome:
    1. Read intervals belonging to this chromosome
    2. Extract non-overlapping regions using vcf_start/vcf_end coordinates
    3. Concatenate appropriate header + filtered variants
    4. Compress with bgzip and index with tabix
    
    Returns: dict mapping chromosome to (snp_vcf_path, indel_vcf_path)
    """
    print(f"\n[Step 3] Merging interval VCFs by chromosome")
    
    check_file(interval_file)
    check_file(snp_header)
    check_file(indel_header)
    check_dir(output_dir)
    
    # Read and group intervals by chromosome
    intervals = read_intervals(interval_file)
    expected_intervals = {f"{chrom}_{reg_start}_{reg_end}" for chrom, reg_start, reg_end, _, _ in intervals}
    validate_vcf_interval_matches(vcf_dir, expected_intervals)
    chr_intervals: Dict[str, List[Tuple[int, int, int, int]]] = {}
    
    for chrom, reg_start, reg_end, vcf_start, vcf_end in intervals:
        if chrom not in chr_intervals:
            chr_intervals[chrom] = []
        chr_intervals[chrom].append((reg_start, reg_end, vcf_start, vcf_end))
    
    print(f"Found {len(intervals)} intervals across {len(chr_intervals)} chromosomes")
    
    output_vcfs = {}
    
    # Process each chromosome
    for chrom in sorted(chr_intervals.keys()):
        print(f"\n[Merging {chrom}]")
        
        # Process SNPs
        vcf_type = 'snp'
        output_vcf = os.path.join(output_dir, f"{chrom}.{vcf_type}.vcf.gz")
        
        if is_file_complete(output_vcf):
            print(f"  Skip merge (exists): {output_vcf}")
            output_vcfs[f"{chrom}_{vcf_type}"] = output_vcf
        else:
            # Build merge command with SNP header
            merge_parts = [f"cat {shlex.quote(snp_header)}"]
            
            for reg_start, reg_end, vcf_start, vcf_end in chr_intervals[chrom]:
                # Construct interval VCF filename
                interval_name = f"{chrom}_{reg_start}_{reg_end}"
                interval_vcf = locate_interval_vcf(vcf_dir, chrom, interval_name, vcf_type)
                if interval_vcf is None:
                    missing_vcf = os.path.join(vcf_dir, chrom, f"{interval_name}.{vcf_type}.vcf.gz")
                    print(f"  Warning: Missing {vcf_type} VCF: {missing_vcf}")
                    continue

                # Add filtered variant extraction
                # Use awk to filter position in vcf_start-vcf_end range (non-overlapping)
                merge_parts.append(f"zcat {shlex.quote(interval_vcf)} | grep -v '^#' | awk '$2>={vcf_start} && $2<={vcf_end}'")
            
            # Combine all parts and compress
            merge_cmd = f"({' ; '.join(merge_parts)}) | {BGZIP} -@ {thread} -c > {output_vcf}"
            
            executor = CommandExecutor(merge_cmd, dry_run, f"Merge {vcf_type.upper()}")
            executor.run()
            
            # Index with tabix
            index_vcf(output_vcf, dry_run)
            
            output_vcfs[f"{chrom}_{vcf_type}"] = output_vcf
        
        # Process INDELs
        vcf_type = 'indel'
        output_vcf = os.path.join(output_dir, f"{chrom}.{vcf_type}.vcf.gz")
        
        if is_file_complete(output_vcf):
            print(f"  Skip merge (exists): {output_vcf}")
            output_vcfs[f"{chrom}_{vcf_type}"] = output_vcf
        else:
            # Build merge command with INDEL header
            merge_parts = [f"cat {shlex.quote(indel_header)}"]
            
            for reg_start, reg_end, vcf_start, vcf_end in chr_intervals[chrom]:
                # Construct interval VCF filename
                interval_name = f"{chrom}_{reg_start}_{reg_end}"
                interval_vcf = locate_interval_vcf(vcf_dir, chrom, interval_name, vcf_type)
                if interval_vcf is None:
                    missing_vcf = os.path.join(vcf_dir, chrom, f"{interval_name}.{vcf_type}.vcf.gz")
                    print(f"  Warning: Missing {vcf_type} VCF: {missing_vcf}")
                    continue

                # Add filtered variant extraction
                # Use awk to filter position in vcf_start-vcf_end range (non-overlapping)
                merge_parts.append(f"zcat {shlex.quote(interval_vcf)} | grep -v '^#' | awk '$2>={vcf_start} && $2<={vcf_end}'")
            
            # Combine all parts and compress
            merge_cmd = f"({' ; '.join(merge_parts)}) | {BGZIP} -@ {thread} -c > {output_vcf}"
            
            executor = CommandExecutor(merge_cmd, dry_run, f"Merge {vcf_type.upper()}")
            executor.run()
            
            # Index with tabix
            index_vcf(output_vcf, dry_run)
            
            output_vcfs[f"{chrom}_{vcf_type}"] = output_vcf
    
    return output_vcfs


# ====================== Main Function ======================
def main():
    parser = argparse.ArgumentParser(
        description='GVCF to VCF Pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Step 0: Generate intervals from fai file
  python %(prog)s step0 -f genome.fa.fai -o intervals.txt

  # Step 1a: Extract GVCF regions (generate command files)
  python %(prog)s step1a -i gvcf_list.txt -l intervals.txt -o interval_gvcfs/ --print-commands

  # Step 1b: Combine GVCFs by interval
  # Mode 1: Process all intervals
  python %(prog)s step1b -i gvcf_list.txt -l intervals.txt -o interval_gvcfs/ -r reference.fa
  
  # Mode 2: Process single interval directory
  python %(prog)s step1b --interval-dir interval_gvcfs/chr1/chr1_1_5000000 -r reference.fa

  # Step 2: Genotype a single GVCF file and generate VCFs
  python %(prog)s step2 -g interval_gvcfs/chr1/chr1_1_5000000.gvcf -o final_vcfs/ -r reference.fa

  # Step 3: Merge interval VCFs by chromosome (using separate headers)
  python %(prog)s step3 -l intervals.txt -v interval_vcfs/ -o final_vcfs/ --snp-header snp_header.txt --indel-header indel_header.txt

  # Step 3: Merge interval VCFs by chromosome (using common header for both)
  python %(prog)s step3 -l intervals.txt -v interval_vcfs/ -o final_vcfs/ -H common_header.txt
        """
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Pipeline steps')
    
    # Step 0: Generate intervals
    step0_parser = subparsers.add_parser('step0', help='Generate intervals from fai file')
    step0_parser.add_argument('-f', '--fai', required=True, help='Genome fai file')
    step0_parser.add_argument('-o', '--output', required=True, help='Output interval file')
    step0_parser.add_argument('-w', '--window', type=int, default=5000000, help='Window size (default: 5000000)')
    step0_parser.add_argument('--overlap', type=int, default=20000, help='Overlap size between windows (default: 20000)')
    
    # Step 1a: Extract GVCF regions
    step1a_parser = subparsers.add_parser('step1a', help='Extract GVCF regions using bcftools')
    step1a_parser.add_argument('-i', '--input', required=True, help='GVCF file list (format: sample_name<TAB>gvcf_path or just gvcf_path)')
    step1a_parser.add_argument('-l', '--intervals', required=True, help='Interval file from step0')
    step1a_parser.add_argument('-o', '--output', required=True, help='Output directory for interval GVCFs')
    step1a_parser.add_argument('-t', '--thread', type=int, default=8, help='Number of threads')
    step1a_parser.add_argument('--dry-run', action='store_true', help='Dry run mode')
    step1a_parser.add_argument('--print-commands', action='store_true', help='Only print commands (one per line) for batch submission')
    
    # Step 1b: Combine GVCFs
    step1b_parser = subparsers.add_parser('step1b', help='Combine extracted GVCFs by interval')
    step1b_parser.add_argument('-i', '--input', help='GVCF file list (required for batch mode)')
    step1b_parser.add_argument('-l', '--intervals', help='Interval file from step0 (required for batch mode)')
    step1b_parser.add_argument('-o', '--output', help='Output directory for interval GVCFs (required for batch mode)')
    step1b_parser.add_argument('--interval-dir', help='Single interval directory to process (alternative to batch mode)')
    step1b_parser.add_argument('-r', '--reference', required=True, help='Reference genome fasta')
    step1b_parser.add_argument('--java-options', type=str, default='-Xmx40g', help='Java options for GATK')
    step1b_parser.add_argument('--dry-run', action='store_true', help='Dry run mode')
    step1b_parser.add_argument('--keep-intermediate', action='store_true', help='Keep intermediate files')
    
    # Step 2: Genotype and generate VCFs
    step2_parser = subparsers.add_parser('step2', help='Genotype a single GVCF file and generate VCFs')
    step2_parser.add_argument('-g', '--gvcf', required=True, help='Input GVCF file')
    step2_parser.add_argument('-o', '--output', required=True, help='Output directory for final VCFs')
    step2_parser.add_argument('-r', '--reference', required=True, help='Reference genome fasta')
    step2_parser.add_argument('-t', '--thread', type=int, default=8, help='Number of threads')
    step2_parser.add_argument('--java-options', type=str, default='-Xmx40g', help='Java options for GATK')
    step2_parser.add_argument('--dry-run', action='store_true', help='Dry run mode')
    step2_parser.add_argument('--keep-intermediate', action='store_true', help='Keep intermediate files')
    step2_parser.add_argument('-p', '--prefix', type=str, default=None, help='Output prefix (default: use input GVCF basename)')
    
    # Step 3: Merge interval VCFs by chromosome
    step3_parser = subparsers.add_parser('step3', help='Merge interval VCFs by chromosome')
    step3_parser.add_argument('-l', '--intervals', required=True, help='Interval file from step0')
    step3_parser.add_argument('-v', '--vcf-dir', required=True, help='Directory containing interval VCFs')
    step3_parser.add_argument('-o', '--output', required=True, help='Output directory for final merged VCFs')
    step3_parser.add_argument('-H', '--header', help='Common header file for both SNP and INDEL (alternative to separate headers)')
    step3_parser.add_argument('--snp-header', help='Header file for SNP VCFs')
    step3_parser.add_argument('--indel-header', help='Header file for INDEL VCFs')
    step3_parser.add_argument('-t', '--thread', type=int, default=8, help='Number of threads for bgzip')
    step3_parser.add_argument('--dry-run', action='store_true', help='Dry run mode')
    
    args = parser.parse_args()
    
    if not args.command:
        parser.print_help()
        sys.exit(1)
    
    # Execute based on command
    if args.command == 'step0':
        generate_intervals(args.fai, args.output, args.window, args.overlap)
    
    elif args.command == 'step1a':
        extract_region_gvcfs(
            args.input,
            args.intervals,
            args.output,
            args.thread,
            args.dry_run,
            args.print_commands if hasattr(args, 'print_commands') else False
        )
    
    elif args.command == 'step1b':
        # Determine mode
        if args.interval_dir:
            # Single interval mode
            combine_region_gvcfs(
                gvcf_list=None,
                interval_file=None,
                output_dir=None,
                reference=args.reference,
                java_options=args.java_options,
                dry_run=args.dry_run,
                keep_intermediate=args.keep_intermediate,
                interval_dir=args.interval_dir
            )
        else:
            # Batch mode - require all parameters
            if not all([args.input, args.intervals, args.output]):
                print("Error: Batch mode requires -i, -l, and -o parameters", file=sys.stderr)
                sys.exit(1)
            combine_region_gvcfs(
                args.input,
                args.intervals,
                args.output,
                args.reference,
                args.java_options,
                args.dry_run,
                args.keep_intermediate
            )
    
    elif args.command == 'step2':
        genotype_and_generate_vcf(
            args.gvcf,
            args.output,
            args.reference,
            args.thread,
            args.java_options,
            args.dry_run,
            args.keep_intermediate,
            args.prefix
        )
    
    elif args.command == 'step3':
        # Validate header parameters
        if not args.header and not (args.snp_header and args.indel_header):
            print("Error: Must provide either --header (common) or both --snp-header and --indel-header", file=sys.stderr)
            sys.exit(1)
        
        # Resolve header files
        snp_header = args.snp_header or args.header
        indel_header = args.indel_header or args.header
        
        merge_interval_vcfs(
            args.intervals,
            args.vcf_dir,
            args.output,
            snp_header,
            indel_header,
            args.thread,
            args.dry_run
        )


if __name__ == '__main__':
    main()