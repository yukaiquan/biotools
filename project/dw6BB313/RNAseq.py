#!/usr/bin/env python
import sys
import os
from subprocess import Popen, PIPE
import time
from concurrent.futures import ProcessPoolExecutor
import click


R_base = '/public/share/acfaa2ssz7/miniconda3/envs/R/bin/Rscript'
Feature_counts = '/public/home/acfaa2ssz7/soft/RNAseqScript/run-featurecounts.R'
Merge_featurecounts = '/public/home/acfaa2ssz7/soft/RNAseqScript/abundance_estimates_to_matrix.pl'
GATK = '/public/home/acfaa2ssz7/soft/RNAseqScript/gatk4.py'

# 判断上方文件是否存在
if not os.path.exists(R_base):
    print('Error: R_base not exist')
    sys.exit(1)
if not os.path.exists(Feature_counts):
    print('Error: Feature_counts not exist')
    sys.exit(1)
if not os.path.exists(Merge_featurecounts):
    print('Error: Merge_featurecounts not exist')
    sys.exit(1)
if not os.path.exists(GATK):
    print('Error: GATK not exist')
    sys.exit(1)


def cacu_threads(threads: int, intheads: int) -> int:
    if threads <= intheads:
        parallel = 1
        intheads = threads
    elif threads > intheads and threads % intheads == 0:
        parallel = threads // intheads
    else:
        parallel = threads // intheads + 1
    return parallel


def run_fastp(fastq1, fastq2, outname1, outname2, threads, outdir):
    # 测试cmd 睡眠10s
    # cmd = 'echo "hello world"\nsleep 2'
    # fastp最大线程数为16
    if outname1 != outname2:
        print('Error: outname1 != outname2')
        sys.exit(1)
    cmd = f'fastp -i {fastq1} -I {fastq2} -o {outdir}/out_{fastq1} -O {outdir}/out_{fastq2} -h {outdir}/{outname1}_fastp.html -j {outdir}/{outname1}_fastp.json -w {threads}'
    p = Popen(cmd, shell=True)
    return_code = p.wait()  # 等待子进程结束，并返回状态码；
    return return_code


def fastqc(fastq1, fastq2, outdir):
    # 测试cmd 睡眠10s
    # cmd = 'echo "hello world"\nsleep 2'
    # fastqc最大线程数设置为8
    cmd = f'fastqc -t 8 -o {outdir} {fastq1} {fastq2}'
    p = Popen(cmd, shell=True)
    return_code = p.wait()  # 等待子进程结束，并返回状态码；
    return return_code


def run_star(genomeDir, fastq1, fastq2, header, outdir, threads):
    # 测试cmd 睡眠10s
    # cmd = 'echo "hello world"\nsleep 2'
    # 如果STAR最大线程数为32 占用最少120G内存
    # 燕麦SFS --alignIntronMin 2 --alignIntronMax 7778 最小内含子长度为2 最大内含子长度为7778
    cmd = f'STAR --twopassMode Basic --quantMode TranscriptomeSAM GeneCounts --runThreadN {threads} --genomeDir {genomeDir} --alignIntronMin 2 --alignIntronMax 7778 --outSAMtype BAM SortedByCoordinate --sjdbOverhang 149 --outSAMattrRGline ID:{header} SM:{header} PL:ILLUMINA --outFilterMismatchNmax 2 --outSJfilterReads Unique --outSAMmultNmax 1 --outFileNamePrefix {outdir} --outSAMmapqUnique 60 --readFilesCommand gunzip -c --readFilesIn {fastq1} {fastq2}'
    p = Popen(cmd, shell=True)
    return_code = p.wait()  # 等待子进程结束，并返回状态码；
    return return_code


def run_featurecounts(R_base, bam, gtf, name, threads):
    cmd = f'{R_base} {Feature_counts} -b {bam} -g {gtf} -o {name} -t {threads}'
    p = Popen(cmd, shell=True)
    return_code = p.wait()  # 等待子进程结束，并返回状态码；
    return return_code


def run_merge_featurecounts(out, sample_list):
    # perl ~/soft/abundance_estimates_to_matrix.pl --est_method featureCounts --cross_sample_norm TMM --out_prefix dw6.matrix --quant_files
    cmd = f'perl {Merge_featurecounts} --est_method featureCounts --cross_sample_norm TMM --out_prefix {out} --quant_files {sample_list}'
    p = Popen(cmd, shell=True)
    return_code = p.wait()  # 等待子进程结束，并返回状态码；
    return return_code


def gatk_remove_duplicates(bam, out):
    cmd = f'gatk MarkDuplicates -I {bam} -O {out} -M {out}.metrics.txt'
    p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    return_code = p.wait()  # 等待子进程结束，并返回状态码；
    stdout, stderr = p.communicate()
    return stdout, stderr


def gatk_add_header(bam, header, out):
    cmd = f'gatk AddOrReplaceReadGroups -I {bam} -O {out} -RGID {header} -RGLB {header} -RGPL ILLUMINA -RGPU unit1 -RGSM {header}'
    p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    return_code = p.wait()  # 等待子进程结束，并返回状态码；
    stdout, stderr = p.communicate()
    return stdout, stderr


def gatk4(ref, input_bam, output):
    cmd = f'python {GATK} -P 21 -I {input_bam} -R {ref} -O {output}'
    p = Popen(cmd, shell=True)
    return_code = p.wait()  # 等待子进程结束，并返回状态码；
    return True


def gatk4_merge(chrom_gvcf, output):
    cmd = f'gatk GatherVcfs -RI TRUE -I {" -I ".join(chrom_gvcf)} -O {output}'
    p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    stdout, stderr = p.communicate()
    return stdout, stderr


def samtools_index(bam):
    cmd = f'samtools index -c {bam}'
    p = Popen(cmd, shell=True)
    return_code = p.wait()
    if return_code != 0:
        print(f'Error: {cmd}')
    return True


def g_vcf_index(gvcf):
    '''
    bcftools version is 1.3.1
    '''
    cmd = f'bcftools index -t {gvcf}'
    if os.path.exists(f'{gvcf}.tbi'):
        return True
    p = Popen(cmd, shell=True)
    return_code = p.wait()
    if return_code != 0:
        print(f'Error: {cmd}')
    return True


def geno_type_GVCFs(ref, input_gvcf, output_vcf):
    cmd = f'gatk GenotypeGVCFs -OVI False -R {ref} -V {input_gvcf} -O {output_vcf}'
    print(cmd)
    if os.path.exists(output_vcf):
        return True
    p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    stdout, stderr = p.communicate()
    return_code = p.wait()
    if return_code != 0:
        print(f'Error: {cmd} {stderr}')
        sys.exit(1)
    return True


def gvcf_merge(ref, gvcf_list, output):
    cmd = f'gatk CombineGVCFs -OVI False -R {ref} -V {" -V ".join(gvcf_list)} -O {output}'
    print(cmd)
    if os.path.exists(output):
        return True
    p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    stdout, stderr = p.communicate()
    return_code = p.wait()
    if return_code != 0:
        print(f'Error: {cmd} {stderr}')
        sys.exit(1)
    return True


def mv_file(infile, outfile):
    cmd = f'mv {infile} {outfile}'
    p = Popen(cmd, shell=True)
    return_code = p.wait()  # 等待子进程结束，并返回状态码；
    return return_code


def mk_sample_dict(sample_list: str) -> dict:
    # 每两行为一组，第一行为R1，第二行为R2
    # 创建数组储存
    sample_dict = {}
    with open(sample_list, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            if line[0].endswith('1.fastq.gz') or line[0].endswith('1.fq.gz'):
                fastq1 = line[0]
                outname1 = line[1]
                NIL = line[2]
                Phenotype = line[3]
            elif line[0].endswith('2.fastq.gz') or line[0].endswith('2.fq.gz'):
                fastq2 = line[0]
                outname2 = line[1]
                NIL = line[2]
                Phenotype = line[3]
                if outname1 == outname2 and NIL == NIL and Phenotype == Phenotype:
                    sample_dict[outname1] = [fastq1, fastq2, NIL, Phenotype]
            else:
                print('Error: fastq name error!')
                sys.exit(1)
    return sample_dict


# 获取参数
@click.command()
@click.option('--sample_list', '-s', help='sample_list')
@click.option('--ref_genome', '-r', help='ref_genome')
@click.option('--star_index', '-i', help='star_index')
@click.option('--gtf', '-g', help='gtf')
@click.option('--threads', '-t', help='threads', default=16)
@click.option('--gatk', '-k', help='gatk', default=False)
def main(sample_list, ref_genome, star_index, gtf, threads, gatk):
    # outdir = '/home/yukaiquan/Ex-seq/data'
    start_time = time.time()
    fastp_outdir = '01_fastp'
    star_outdir = '02_star'
    feature_outdir = '03_featureCounts'
    matrix_outdir = '04_matrix'
    add_outdir = '05_add'
    dup_bam_outdir = '06_dup_bam'
    gvcf_outdir = '07_gvcf'
    vcf_outdir = '08_vcf'
    # 创建输出目录
    if not os.path.exists(fastp_outdir):
        os.makedirs(fastp_outdir)
    sample_dict = mk_sample_dict(sample_list)
    # run_fastp(fastq1, fastq2, outdir, threads)
    # 整理sample_list
    fastp_theads = 16
    parallel = cacu_threads(threads, fastp_theads)
    with ProcessPoolExecutor(max_workers=parallel) as executor:
        for outname1, [fastq1, fastq2, NIL, _] in sample_dict.items():
            # print(outname1, fastq1, fastq2)
            # executor.submit(run_fastp, fastq1, fastq2, outname1,
            #                 outname1, fastp_theads, fastp_outdir)
            print(f'{outname1} fastp done!')
    print('#'*25, 'fastp done!', '#'*25)
    # 判断线程数计算并行数
    # fastq_theads = 8
    # parallel = cacu_threads(threads, fastq_theads)
    # with ProcessPoolExecutor(max_workers=parallel) as executor:
    #     for outname1, [fastq1, fastq2, NIL, _] in sample_dict.items():
    #         executor.submit(fastqc, f'{fastp_outdir}/out_{fastq1}',
    #                         f'{fastp_outdir}/out_{fastq2}', fastp_outdir)
    # print('#'*25, 'fastqc done!', '#'*25)
    if not os.path.exists(star_outdir):
        os.makedirs(star_outdir)
    # STAR在不满节点时,如果开并行,会导致内存不足,所以不开并行
    for outname1, [fastq1, fastq2, NIL, _] in sample_dict.items():
        if run_star(f'{star_index}', f'{fastp_outdir}/out_{fastq1}', f'{fastp_outdir}/out_{fastq2}', f'{NIL}', f'{star_outdir}' + '/' + f'{outname1}', f'{threads}') == 0:
            print(f'{outname1} STAR done!')
        else:
            print(f'{outname1} STAR error!')
            sys.exit(1)
    print('#'*25, 'STAR done!', '#'*25)
    if not os.path.exists(feature_outdir):
        os.makedirs(feature_outdir)
    # 判断线程数计算并行数
    feature_theads = 8
    parallel = cacu_threads(threads, feature_theads)
    with ProcessPoolExecutor(max_workers=parallel) as executor:
        for outname1, [fastq1, fastq2, NIL, Phenotype] in sample_dict.items():
            executor.submit(
                run_featurecounts, f'{R_base}', f'{star_outdir}/{outname1}Aligned.sortedByCoord.out.bam', f'{gtf}', f'{feature_outdir}/{outname1}', f'{feature_theads}')
            print(f'{outname1} featurecounts done!')

    # 更改输出的文件名字# 提取NIL相同的样本
    NIL_dict = {}
    for outname1, [fastq1, fastq2, NIL, Phenotype] in sample_dict.items():
        if os.path.exists(f'{feature_outdir}/{outname1}.count') and os.path.exists(f'{feature_outdir}/{outname1}.log'):
            mv_file(f'{feature_outdir}/{outname1}.count',
                    f'{feature_outdir}/{outname1}_{NIL}_{Phenotype}.count')
            mv_file(f'{feature_outdir}/{outname1}.log',
                    f'{feature_outdir}/{outname1}_{NIL}_{Phenotype}.log')
        if NIL not in NIL_dict:
            NIL_dict[NIL] = [f'{outname1}_{NIL}_{Phenotype}']
        else:
            NIL_dict[NIL].append(f'{outname1}_{NIL}_{Phenotype}')

    # 根据NIL dict写出文件
    for NIL, sample_list in NIL_dict.items():
        with open(f'{feature_outdir}/{NIL}.txt', 'w') as f:
            for sample in sample_list:
                f.write(f'{feature_outdir}/{sample}.count\n')
    # 新建目录保存结果
    if not os.path.exists(matrix_outdir):
        os.makedirs(matrix_outdir)

    for NIL, sample_list in NIL_dict.items():
        if run_merge_featurecounts(f'{NIL}', f'{feature_outdir}/{NIL}.txt') == 0:
            # 移动结果保存
            mv_file(f'{NIL}.*', f'{matrix_outdir}/')
            print(f'{NIL} merge Done !')
        else:
            print('Error')
            sys.exit(1)

    if gatk:
        #######################################call SNP############################################
        if not os.path.exists(add_outdir):
            os.makedirs(add_outdir)
        if not os.path.exists(dup_bam_outdir):
            os.makedirs(dup_bam_outdir)
        # 判断线程数计算并行数
        dup_theads = 8
        parallel = cacu_threads(threads, dup_theads)
        with ProcessPoolExecutor(max_workers=parallel) as executor:
            for outname1, [fastq1, fastq2, NIL, Phenotype] in sample_dict.items():
                executor.submit(
                    gatk_add_header, f'{star_outdir}/{outname1}Aligned.sortedByCoord.out.bam', f'{NIL}{Phenotype}', f'{add_outdir}/{outname1}_header.bam')
        # 判断线程数计算并行数
        add_theads = 8
        parallel = cacu_threads(threads, add_theads)
        with ProcessPoolExecutor(max_workers=parallel) as executor:
            for outname1, [fastq1, fastq2, NIL, Phenotype] in sample_dict.items():
                executor.submit(
                    gatk_remove_duplicates, f'{add_outdir}/{outname1}_header.bam', f'{dup_bam_outdir}/{outname1}_dup_add.bam')
        # 判断线程数计算并行数
        samtools_theads = 2
        parallel = cacu_threads(threads, samtools_theads)
        # 创建文件夹
        if not os.path.exists(gvcf_outdir):
            os.makedirs(gvcf_outdir)
        with ProcessPoolExecutor(max_workers=parallel) as executor:
            for outname1, [fastq1, fastq2, NIL, Phenotype] in sample_dict.items():
                if not os.path.exists(f'{dup_bam_outdir}/{outname1}_dup_add.bam.csi'):
                    executor.submit(
                        samtools_index, f'{dup_bam_outdir}/{outname1}_dup_add.bam')
                print(f'{outname1} gvcf to vcf Done !')
        # gatk 需要21线程 且需要索引CSI 最好一次在一个节点调用一次需要提供超过21线程的计算资源
        try:
            for outname1, [fastq1, fastq2, NIL, Phenotype] in sample_dict.items():
                if not os.path.exists(f'{gvcf_outdir}/{outname1}.g.vcf.gz'):
                    gatk4(f'{ref_genome}',
                          f'{dup_bam_outdir}/{outname1}_dup_add.bam', f'{outname1}')
        except Exception as e:
            print(e)
        # 移动文件 建立索引
        for outname1, [fastq1, fastq2, NIL, Phenotype] in sample_dict.items():
            if not os.path.exists(f'{gvcf_outdir}/{outname1}.g.vcf.gz'):
                mv_file(f'{outname1}.*g.vcf.gz', f'{gvcf_outdir}/')
            if not os.path.exists(f'{gvcf_outdir}/{outname1}.g.vcf.gz.tbi'):
                g_vcf_index(f'{gvcf_outdir}/{outname1}.g.vcf.gz')
        # 产生vcf文件存储位置
        if not os.path.exists(vcf_outdir):
            os.makedirs(vcf_outdir)
        # 合并gvcf
        gvcf_dict = {}
        for outname1, [fastq1, fastq2, NIL, Phenotype] in sample_dict.items():
            # 所有样品合并gvcfcall snp
            # if 'All_samples' not in gvcf_dict:
            #     gvcf_dict['All_samples'] = [
            #         f'{gvcf_outdir}/{outname1}.g.vcf.gz']
            # else:
            #     gvcf_dict['All_samples'].append(
            #         f'{gvcf_outdir}/{outname1}.g.vcf.gz')
            if f'{NIL}{Phenotype}' not in gvcf_dict:
                gvcf_dict[f'{NIL}{Phenotype}'] = [
                    f'{gvcf_outdir}/{outname1}.g.vcf.gz']
            else:
                gvcf_dict[f'{NIL}{Phenotype}'].append(
                    f'{gvcf_outdir}/{outname1}.g.vcf.gz')
        # 对字典中列表进行去重 防止重复运行
        gvcf_dict = {k: list(set(v)) for k, v in gvcf_dict.items()}
        parallel = cacu_threads(threads, 4)
        with ProcessPoolExecutor(max_workers=parallel) as executor:
            for header, gvcf_list in gvcf_dict.items():
                if not os.path.exists(f'{vcf_outdir}/{header}.g.vcf.gz'):
                    executor.submit(
                        gvcf_merge, ref_genome, gvcf_list, f'{vcf_outdir}/{header}.g.vcf.gz')
        # 计算并行数
        parallel = cacu_threads(threads, 1)
        with ProcessPoolExecutor(max_workers=parallel) as executor:
            for header, gvcf_list in gvcf_dict.items():
                if not os.path.exists(f'{vcf_outdir}/{header}.g.vcf.gz.tbi'):
                    executor.submit(
                        g_vcf_index, f'{vcf_outdir}/{header}.g.vcf.gz')
        # 计算并行数
        parallel = cacu_threads(threads, 4)
        with ProcessPoolExecutor(max_workers=parallel) as executor:
            for header, gvcf_list in gvcf_dict.items():
                if not os.path.exists(f'{vcf_outdir}/{header}.vcf.gz'):
                    executor.submit(
                        geno_type_GVCFs, ref_genome, f'{vcf_outdir}/{header}.g.vcf.gz', f'{vcf_outdir}/{header}.vcf.gz')
        # 计算并行数
        parallel = cacu_threads(threads, 1)
        with ProcessPoolExecutor(max_workers=parallel) as executor:
            for header, gvcf_list in gvcf_dict.items():
                if not os.path.exists(f'{vcf_outdir}/{header}.vcf.gz.tbi'):
                    executor.submit(
                        g_vcf_index, f'{vcf_outdir}/{header}.vcf.gz')
    ###########################################################
    print('All Done !')
    end_time = time.time()
    print('RNAseq Time used: %s minute' % (end_time - start_time / 60))


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print('KeyboardInterrupt')
        sys.exit(0)
