#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
import os
import subprocess
import sys
import glob
import datetime  # 导入时间模块，用于日志时间戳


# SIF = "singularity exec ~/soft/kmergwas.sif"
KMC = f"kmc"
KMC_STRAND = f"/public/share/h13713/soft/kmersGWAS/bin/kmers_add_strand_information"


def get_fastq_files(folder):
    """获取文件夹中的所有双端 fastq.gz 文件路径"""
    fastq_files = []
    for root, _, files in os.walk(folder):
        for file in files:
            # 筛选双端文件（含_1/_2标记）和fastq.gz/fq.gz后缀
            if (file.endswith("fq.gz") or file.endswith("fastq.gz")) and ("1" in file or "2" in file):
                fastq_files.append(os.path.join(root, file))
    return fastq_files


def write_file_list(file_list, output_path, log_path):
    """将 fastq 文件路径写入 file.list，并记录日志"""
    with open(output_path, "w") as f:
        for file in file_list:
            f.write(f"{file}\n")
    # 记录日志：生成file.list的结果
    current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    log_entry = f"[{current_time}] 已生成fastq文件列表：{output_path}\n" \
                f"[{current_time}] 共包含 {len(file_list)} 个双端fastq文件\n"
    print(log_entry.strip())  # 控制台同步输出
    with open(log_path, "a", encoding="utf-8") as f:
        f.write(log_entry)


def run_command(command, log_path):
    """运行命令并记录详细日志（时间戳、命令内容、stdout/stderr、执行结果）"""
    current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    # 记录“开始执行命令”的日志
    start_log = f"[{current_time}] 开始执行命令：\n{command}\n"
    print(f"[{current_time}] 开始执行命令...")  # 控制台简化提示
    with open(log_path, "a", encoding="utf-8") as f:
        f.write(start_log)
    
    try:
        # 捕获命令的stdout和stderr，方便日志记录
        result = subprocess.run(
            command, 
            shell=True, 
            check=True, 
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE,
            universal_newlines=True  # 替换为Python 3.6支持的等效参数
        )
        # 记录“命令执行成功”的日志（含stdout）
        success_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        success_log = f"[{success_time}] 命令执行成功！\n" \
                     f"[{success_time}] 命令stdout输出：\n{result.stdout}\n"
        print(f"[{success_time}] 命令执行成功")  # 控制台简化提示
        with open(log_path, "a", encoding="utf-8") as f:
            f.write(success_log)
    
    except subprocess.CalledProcessError as e:
        # 记录“命令执行失败”的日志（含stdout/stderr和错误信息）
        error_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        error_log = f"[{error_time}] 命令执行失败！\n" \
                   f"[{error_time}] 错误信息：{str(e)}\n" \
                   f"[{error_time}] 命令stdout输出：\n{e.stdout}\n" \
                   f"[{error_time}] 命令stderr输出：\n{e.stderr}\n"
        print(error_log.strip())  # 控制台打印完整错误
        with open(log_path, "a", encoding="utf-8") as f:
            f.write(error_log)
        raise  # 重新抛出异常，触发脚本退出（确保错误不被忽略）


def delete_files_with_prefix(folder, prefix, log_path):
    """删除文件夹中以指定前缀开头的文件，并记录删除日志"""
    files_to_delete = glob.glob(os.path.join(folder, f"{prefix}*"))
    current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    if not files_to_delete:
        log_entry = f"[{current_time}] 未找到以「{prefix}」为前缀的文件，无需删除\n"
        print(log_entry.strip())
        with open(log_path, "a", encoding="utf-8") as f:
            f.write(log_entry)
        return
    
    # 逐个删除文件并记录日志
    for file in files_to_delete:
        try:
            os.remove(file)
            delete_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            delete_log = f"[{delete_time}] 已删除文件：{file}\n"
            print(delete_log.strip())
            with open(log_path, "a", encoding="utf-8") as f:
                f.write(delete_log)
        except OSError as e:
            error_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            error_log = f"[{error_time}] 删除文件失败：{file}\n" \
                       f"[{error_time}] 错误信息：{str(e)}\n"
            print(error_log.strip())
            with open(log_path, "a", encoding="utf-8") as f:
                f.write(error_log)


def check_step_output(output_files, step_name, log_path):
    """
    检查步骤输出文件是否存在且有效
    :param output_files: 输出文件路径列表（可以是单个文件或多个文件）
    :param step_name: 步骤名称（用于日志）
    :param log_path: 日志文件路径
    :return: True（文件存在且有效，跳过步骤）/ False（文件不存在或无效，需要执行步骤）
    """
    current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    # 确保output_files是列表
    if not isinstance(output_files, list):
        output_files = [output_files]
    
    # 检查所有文件是否存在且大小>10字节（避免空文件）
    all_valid = True
    for file in output_files:
        if not os.path.exists(file):
            all_valid = False
            log_entry = f"[{current_time}] 步骤「{step_name}」：未找到输出文件 {file}\n"
            print(log_entry.strip())
            with open(log_path, "a", encoding="utf-8") as f:
                f.write(log_entry)
        elif os.path.getsize(file) <= 10:
            all_valid = False
            log_entry = f"[{current_time}] 步骤「{step_name}」：输出文件 {file} 为空（大小≤10字节）\n"
            print(log_entry.strip())
            with open(log_path, "a", encoding="utf-8") as f:
                f.write(log_entry)
    
    if all_valid:
        log_entry = f"[{current_time}] 步骤「{step_name}」：所有输出文件均存在且有效，跳过该步骤\n"
        print(log_entry.strip())
        with open(log_path, "a", encoding="utf-8") as f:
            f.write(log_entry)
        return True
    else:
        log_entry = f"[{current_time}] 步骤「{step_name}」：输出文件不完整或无效，将执行该步骤\n"
        print(log_entry.strip())
        with open(log_path, "a", encoding="utf-8") as f:
            f.write(log_entry)
        return False


def main(folder, threads):
    # 1. 初始化日志文件（路径：目标文件夹/kmc_run.log）
    log_path = os.path.join(folder, "kmc_run.log")
    init_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    # 日志头部：记录脚本启动信息
    with open(log_path, "w", encoding="utf-8") as f:
        f.write(f"[{init_time}] KMC流程脚本启动\n")
        f.write(f"[{init_time}] 目标文件夹：{folder}\n")
        f.write(f"[{init_time}] 使用线程数：{threads}\n\n")
    print(f"[{init_time}] KMC流程启动，日志文件已保存至：{log_path}")

    # 2. 获取fastq文件并生成file.list（步骤1：生成文件列表）
    fastq_files = get_fastq_files(folder)
    if not fastq_files:
        error_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        error_msg = f"[{error_time}] 错误：未在文件夹「{folder}」中找到双端fastq.gz/fq.gz文件\n"
        print(error_msg.strip())
        with open(log_path, "a", encoding="utf-8") as f:
            f.write(error_msg)
        sys.exit(1)
    
    file_list_path = os.path.join(folder, "file.list")
    # 检查file.list是否存在且有效
    if not check_step_output(file_list_path, "生成fastq文件列表", log_path):
        write_file_list(fastq_files, file_list_path, log_path)

    # 3. 定义核心参数（与原逻辑一致）
    kmer_size = 31
    thread_count = threads
    canon_kmer_cutoff = 2  # 规范化k-mer计数阈值（至少出现2次）
    all_kmer_cutoff = 0    # 非规范化k-mer计数阈值（全部统计）
    kmc_canon_output = os.path.join(folder, "output_kmc_canon")
    kmc_all_output = os.path.join(folder, "output_kmc_all")
    kmers_with_strand_output = os.path.join(folder, "kmers_with_strand")

    # KMC命令的输出文件（KMC会生成.kmc_pre和.kmc_suf两个文件）
    kmc_canon_files = [f"{kmc_canon_output}.kmc_pre", f"{kmc_canon_output}.kmc_suf"]
    kmc_all_files = [f"{kmc_all_output}.kmc_pre", f"{kmc_all_output}.kmc_suf"]

    # 4. 构造并执行KMC相关命令（带日志记录和跳过逻辑）
    current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    log_entry = f"\n[{current_time}] 开始执行KMC相关命令（k-mer长度：{kmer_size}）\n"
    print(log_entry.strip())
    with open(log_path, "a", encoding="utf-8") as f:
        f.write(log_entry)

    # 4.1 执行规范化KMC命令（-ci{canon_kmer_cutoff}）
    kmc_command_canon = f"{KMC} -t{thread_count} -k{kmer_size} -ci{canon_kmer_cutoff} @{file_list_path} {kmc_canon_output} {folder}"
    if not check_step_output(kmc_canon_files, "规范化KMC计数", log_path):
        run_command(kmc_command_canon, log_path)

    # 4.2 执行非规范化KMC命令（-b -ci{all_kmer_cutoff}）
    kmc_command_all = f"{KMC} -t{thread_count} -k{kmer_size} -ci{all_kmer_cutoff} -b @{file_list_path} {kmc_all_output} {folder}"
    if not check_step_output(kmc_all_files, "非规范化KMC计数", log_path):
        run_command(kmc_command_all, log_path)

    # 4.3 执行k-mer链信息合并命令
    kmers_add_strand_command = f"{KMC_STRAND} -c {kmc_canon_output} -n {kmc_all_output} -k {kmer_size} -o {kmers_with_strand_output}"
    if not check_step_output(kmers_with_strand_output, "k-mer链信息合并", log_path):
        run_command(kmers_add_strand_command, log_path)

    # 5. 检查输出文件有效性（带日志记录）
    current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    log_entry = f"\n[{current_time}] 开始检查最终输出文件有效性\n"
    print(log_entry.strip())
    with open(log_path, "a", encoding="utf-8") as f:
        f.write(log_entry)

    if os.path.exists(kmers_with_strand_output) and os.path.getsize(kmers_with_strand_output) > 10:
        # 输出文件正常（非空）
        ok_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        ok_log = f"[{ok_time}] 最终输出文件正常：{kmers_with_strand_output}\n" \
                 f"[{ok_time}] 文件大小：{os.path.getsize(kmers_with_strand_output)} 字节\n"
        print(ok_log.strip())
        with open(log_path, "a", encoding="utf-8") as f:
            f.write(ok_log)

        # （可选）删除KMC中间文件（如需启用，取消注释下方两行）
        delete_files_with_prefix(folder, "output_kmc_all", log_path)
        delete_files_with_prefix(folder, "output_kmc_canon", log_path)

        # （可选）删除fastq文件（如需启用，取消注释下方代码块）
        # del_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        # log_entry = f"[{del_time}] 开始删除原始fastq文件\n"
        # print(log_entry.strip())
        # with open(log_path, "a", encoding="utf-8") as f:
        #     f.write(log_entry)
        # for file in fastq_files:
        #     try:
        #         os.remove(file)
        #         del_log = f"[{del_time}] 已删除fastq文件：{file}\n"
        #         print(del_log.strip())
        #         with open(log_path, "a", encoding="utf-8") as f:
        #             f.write(del_log)
        #     except OSError as e:
        #         del_err_log = f"[{del_time}] 删除fastq文件失败：{file}，错误：{str(e)}\n"
        #         print(del_err_log.strip())
        #         with open(log_path, "a", encoding="utf-8") as f:
        #             f.write(del_err_log)

    else:
        # 输出文件为空或不存在（报错并退出）
        err_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        err_log = f"[{err_time}] 错误：最终输出文件异常！\n" \
                 f"[{err_time}] 文件路径：{kmers_with_strand_output}\n" \
                 f"[{err_time}] 异常原因：文件不存在 或 文件大小≤10字节（为空）\n"
        print(err_log.strip())
        with open(log_path, "a", encoding="utf-8") as f:
            f.write(err_log)
        sys.exit(1)

    # 6. 流程结束日志
    end_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    end_log = f"\n[{end_time}] 所有命令执行完成！流程正常结束\n" \
              f"[{end_time}] 完整日志已保存至：{log_path}\n"
    print(end_log.strip())
    with open(log_path, "a", encoding="utf-8") as f:
        f.write(end_log)


if __name__ == "__main__":
    # 命令行参数校验（需传入「文件夹路径」和「线程数」）
    if len(sys.argv) != 3:
        print("用法: python run_kmc.py <文件夹路径> <线程数>")
        print("示例: python run_kmc.py ./sample1 8")
        sys.exit(1)

    folder = sys.argv[1]
    threads = sys.argv[2]

    # 校验文件夹是否存在
    if not os.path.exists(folder):
        print(f"错误：文件夹「{folder}」不存在，请检查路径是否正确")
        sys.exit(1)

    # 启动主流程
    main(folder, threads)
