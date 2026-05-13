import pysam
import argparse
import os

def extract_snps_from_vcf(vcf_input, snp_list_file, vcf_output):
    # 校验输入文件是否存在
    if not os.path.exists(vcf_input):
        raise FileNotFoundError(f"错误：输入VCF文件不存在 → {vcf_input}")
    if not os.path.exists(snp_list_file):
        raise FileNotFoundError(f"错误：SNP列表文件不存在 → {snp_list_file}")
    
    # 1. 读取SNP ID列表，转为集合（查询更快）
    snp_set = set()
    with open(snp_list_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            snp_id = line.strip()
            if snp_id:  # 跳过空行
                snp_set.add(snp_id)
            else:
                print(f"警告：SNP列表第{line_num}行为空，已跳过")
    
    if not snp_set:
        raise ValueError("错误：SNP列表文件中未找到有效SNP ID（可能全为空行）")
    
    print(f"成功读取SNP列表：共 {len(snp_set)} 个有效SNP ID")
    
    # 2. 打开输入VCF和输出VCF（支持.vcf和.vcf.gz）
    try:
        with pysam.VariantFile(vcf_input, 'r') as vcf_in, \
             pysam.VariantFile(vcf_output, 'w', header=vcf_in.header) as vcf_out:
            
            # 3. 遍历VCF，提取目标SNP
            extracted_count = 0
            for record in vcf_in:
                # 处理ID字段：支持单个ID、多个ID（;分隔）、无ID（.）
                record_ids = record.id.split(';') if (record.id and record.id != '.') else []
                if any(snp in snp_set for snp in record_ids):
                    vcf_out.write(record)
                    extracted_count += 1
        
        # 输出统计信息
        print("\n" + "="*50)
        print(f"提取完成！结果文件：{vcf_output}")
        print(f"目标SNP总数：{len(snp_set)}")
        print(f"成功提取：{extracted_count} 个")
        print(f"未找到的SNP：{len(snp_set) - extracted_count} 个")
        print("="*50)
    
    except Exception as e:
        raise RuntimeError(f"处理VCF文件时出错：{str(e)}")

def main():
    # 创建命令行参数解析器
    parser = argparse.ArgumentParser(
        description="从VCF文件中提取指定SNP ID列表对应的SNP集合（支持大文件和压缩格式）",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="""
使用示例：
  1. 基本用法（未压缩VCF）：
     python extract_snps_cli.py -i input.vcf -s snp_list.txt -o extracted.vcf
  
  2. 处理压缩VCF（.vcf.gz）：
     python extract_snps_cli.py -i input.vcf.gz -s snp_list.txt -o extracted.vcf.gz
  
  3. 完整路径示例：
     python extract_snps_cli.py \
       -i /data/raw/genome.vcf.gz \
       -s /data/list/target_snps.txt \
       -o /data/result/filtered_snps.vcf
        """
    )
    
    # 添加必需参数
    parser.add_argument(
        "-i", "--input", 
        required=True, 
        help="输入VCF文件路径（支持 .vcf 或 .vcf.gz 格式）"
    )
    parser.add_argument(
        "-s", "--snp-list", 
        required=True, 
        help="SNP ID列表文件路径（每行1个ID，支持空行自动跳过）"
    )
    parser.add_argument(
        "-o", "--output", 
        required=True, 
        help="输出VCF文件路径（可保存为 .vcf 或 .vcf.gz，自动兼容格式）"
    )
    
    # 解析参数
    args = parser.parse_args()
    
    # 执行提取
    try:
        extract_snps_from_vcf(args.input, args.snp_list, args.output)
    except Exception as e:
        print(f"\n❌ 错误：{str(e)}")
        exit(1)

if __name__ == "__main__":
    main()
