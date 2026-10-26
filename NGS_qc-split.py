import os
import gzip
import numpy as np
from Bio import SeqIO
import argparse

# 质量阈值，用于判断低质量碱基
QUALITY_THRESHOLD = 30

def load_barcode_list(barcode_file):
    """加载barcode列表，返回一个包含所有barcode对的集合"""
    with open(barcode_file, 'r') as f:
        barcodes = [line.strip().split() for line in f.readlines()]
    return set(tuple(barcode) for barcode in barcodes)

def is_high_quality(sequence, quality_list, threshold=QUALITY_THRESHOLD):
    """检查一个序列在指定区域内的质量是否高于给定阈值"""
    # 直接使用质量分数列表，无需使用 ord()
    quality_scores = np.array(quality_list)  # Phred质量分数
    return np.all(quality_scores >= threshold)

def check_sgRNA_quality(r1_seq, r1_qual, r2_seq, r2_qual, r1_sgRNA_start, r2_sgRNA_start, sgRNA_length=20, threshold=QUALITY_THRESHOLD):
    """检查 reads 中 sgRNA 区域的质量"""
    # 提取正向和反向 reads 的 sgRNA 区域
    r1_sgRNA_seq = r1_seq[r1_sgRNA_start:r1_sgRNA_start+sgRNA_length]
    r1_sgRNA_qual = r1_qual[r1_sgRNA_start:r1_sgRNA_start+sgRNA_length]

    r2_sgRNA_seq = r2_seq[r2_sgRNA_start:r2_sgRNA_start+sgRNA_length]
    r2_sgRNA_qual = r2_qual[r2_sgRNA_start:r2_sgRNA_start+sgRNA_length]

    # 检查 sgRNA 区域的质量
    r1_sgRNA_quality_ok = is_high_quality(r1_sgRNA_seq, r1_sgRNA_qual, threshold)
    r2_sgRNA_quality_ok = is_high_quality(r2_sgRNA_seq, r2_sgRNA_qual, threshold)

    return r1_sgRNA_quality_ok and r2_sgRNA_quality_ok

def process_reads(r1_record, r2_record, barcode_set, r1_sgRNA_start, r2_sgRNA_start, sgRNA_length=20, threshold=QUALITY_THRESHOLD):
    """处理正向和反向 read，进行质量检查并过滤不合格的 reads"""

    # 提取序列和质量信息
    r1_seq, r1_qual = str(r1_record.seq), r1_record.letter_annotations["phred_quality"]
    r2_seq, r2_qual = str(r2_record.seq), r2_record.letter_annotations["phred_quality"]

    # 第一步：检查 sgRNA 区域的质量
    if not check_sgRNA_quality(r1_seq, r1_qual, r2_seq, r2_qual, r1_sgRNA_start, r2_sgRNA_start, sgRNA_length, threshold):
        return None

    # 第二步：检查前三个碱基的质量
    if not (is_high_quality(r1_seq[:3], r1_qual[:3], threshold) and
            is_high_quality(r2_seq[:3], r2_qual[:3], threshold)):
        return None

    # 第三步：检查 barcode 匹配
    barcode_r1 = r1_seq[:3]
    barcode_r2 = r2_seq[:3]

    # 将barcode_r1和barcode_r2组合成一个元组
    barcode_pair = (barcode_r1, barcode_r2)

    # 检查组合的barcode是否在barcode_set中
    if barcode_pair not in barcode_set:
        return None  # 如果barcode不匹配，丢弃当前这一对reads

    # 如果通过了所有条件，返回处理过的 reads
    return (r1_record, r2_record)

def main(r1_file, r2_file, barcode_file, output_dir):
    """主函数，加载数据并处理"""
    # 加载barcode列表
    barcode_set = load_barcode_list(barcode_file)

    # 创建输出目录
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # 字典保存已打开的文件对象
    output_files = {}

    # 打开输入文件并读取
    with gzip.open(r1_file, 'rt') as r1_f, gzip.open(r2_file, 'rt') as r2_f:
        r1_records = SeqIO.parse(r1_f, "fastq")
        r2_records = SeqIO.parse(r2_f, "fastq")

        # 假设每一对read是按顺序匹配的
        for r1_record, r2_record in zip(r1_records, r2_records):
            # 这里提供sgRNA的起始位置，假设是85和103
            r1_sgRNA_start = 84
            r2_sgRNA_start = 102

            # 处理每对reads
            result = process_reads(r1_record, r2_record, barcode_set, r1_sgRNA_start, r2_sgRNA_start)

            # 如果结果有效（不是None），则拆分到相应的文件中
            if result:
                barcode_r1 = str(result[0].seq[:3])  # 正向barcode
                barcode_r2 = str(result[1].seq[:3])  # 反向barcode
                barcode_pair = (barcode_r1, barcode_r2)

                # 为当前barcode_pair生成文件名，保持.gz格式
                r1_output_filename = f"{os.path.basename(r1_file).split('.')[0]}_QC_{barcode_r1}_{barcode_r2}.fq.gz"
                r2_output_filename = f"{os.path.basename(r2_file).split('.')[0]}_QC_{barcode_r1}_{barcode_r2}.fq.gz"

                # 生成正向和反向文件的完整路径
                r1_output_filepath = os.path.join(output_dir, r1_output_filename)
                r2_output_filepath = os.path.join(output_dir, r2_output_filename)

                # 检查该barcode_pair的输出文件是否已经创建
                if barcode_pair not in output_files:
                    # 如果未创建，创建对应的文件并打开
                    output_files[barcode_pair] = {
                        'r1': gzip.open(r1_output_filepath, 'wt'),
                        'r2': gzip.open(r2_output_filepath, 'wt')
                    }

                # 写入正向和反向reads到对应的gzip文件
                SeqIO.write(result[0], output_files[barcode_pair]['r1'], "fastq")
                SeqIO.write(result[1], output_files[barcode_pair]['r2'], "fastq")

    # 关闭所有输出文件
    for f in output_files.values():
        f['r1'].close()
        f['r2'].close()

# 通过argparse接收命令行参数
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process sequencing reads")
    parser.add_argument('r1_file', help="Path to the R1 (forward) sequencing file")
    parser.add_argument('r2_file', help="Path to the R2 (reverse) sequencing file")
    parser.add_argument('barcode_file', help="Path to the barcode list file")
    parser.add_argument('output_dir', help="Directory to save the output files")

    # 解析命令行参数
    args = parser.parse_args()

    # 执行主函数，传入命令行参数
    main(args.r1_file, args.r2_file, args.barcode_file, args.output_dir)
