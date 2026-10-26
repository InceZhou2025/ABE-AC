# -*- coding: utf-8 -*-

import os
import shutil
import subprocess

# 固定参数
barcode_file = "/NAS/home/yfwang/NGS-DATA/TDN3-ALL/TDN3-10-ALL_qc-split/barcodes.txt"
input_dir = "/NAS/home/yfwang/NGS-DATA/TDN3-ALL/TDN3-10-ALL_qc-split"
output_root = "/NAS/home/yfwang/NGS-DATA/TDN3-ALL/TDN3-10-ALL_CRISPResso2_Output"
amplicon_seq = "tgaatggattccttggaaacaatgataacaagacctggctgagctaactgtgacagcatgtggtaattttccagcccgctggccctgtaaaggaaactggaacacaaagcatagactgcggggcgggccagcctgaatagctgcaaacaagtgcagaatatctgatgatgtcatacgcacagtttgacagatggggctgg"
guide_seq = "gaacacaaagcatagactgc"

# 四个需要保存的输出文件
categories = {
    "BE_percentage": "Quantification_window_nucleotide_percentage_table.txt",
    "BE_frequency": "Quantification_window_nucleotide_frequency_table.txt",
    "PDF": "9.Alleles_frequency_table_around_sgRNA_gaacacaaagcatagactgc.pdf",
    "Indel": "CRISPResso_quantification_of_editing_frequency.txt"
}

# 创建输出子目录
for folder in categories:
    os.makedirs(os.path.join(output_root, folder), exist_ok=True)

# 读取 barcodes
with open(barcode_file) as f:
    barcodes = [line.strip().replace(" ", "_") for line in f if line.strip()]

for idx, barcode in enumerate(barcodes, 1):
    sample_name = f"TDN3-10-ALL_{barcode}"
    r1 = os.path.join(input_dir, f"TDN3-10-ALL_1_QC_{barcode}.fq.gz")
    r2 = os.path.join(input_dir, f"TDN3-10-ALL_2_QC_{barcode}.fq.gz")
    out_dir = os.path.join(output_root, barcode)
    os.makedirs(out_dir, exist_ok=True)

    # 运行CRISPResso
    subprocess.run([
        "CRISPResso",
        "--fastq_r1", r1,
        "--fastq_r2", r2,
        "--amplicon_seq", amplicon_seq,
        "--guide_seq", guide_seq,
        "--output_folder", out_dir,
        "--name", sample_name,
        "--plot_window_size", "30",
        "--quantification_window_size", "10",
        "--quantification_window_center", "-10",
        "--min_average_read_quality", "20",
        "--base_editor_output"
    ], check=True)

    result_dir = os.path.join(out_dir, f"CRISPResso_on_{sample_name}")

    # 移动并重命名目标输出文件
    for cat, fname in categories.items():
        src = os.path.join(result_dir, fname)
        if os.path.isfile(src):
            ext = os.path.splitext(fname)[1]  # 提取原始文件的扩展名，例如 ".pdf"、".txt"
            new_name = f"{idx}_{sample_name}_{cat}{ext}"  # 添加扩展名
            dst = os.path.join(output_root, cat, new_name)
            shutil.copy(src, dst)

    # 删除整个 CRISPResso 输出目录（如 GAG_ACC）
    shutil.rmtree(os.path.dirname(result_dir))
