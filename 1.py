import os
import sys

# 定义输入和输出目录
dir = "./Input_data/"
outdir = "./ATAC_peaks/"
os.makedirs(outdir, exist_ok=True)

# 读取命令行参数
c = sys.argv[1]
reference = sys.argv[2]

print(c)

# 处理输入CSV文件并生成两个临时文件
with open(f"{dir}{c}.csv") as f1, \
     open(f"{outdir}{c}.clean.ATAC_peak.temp", "w") as o1, \
     open(f"{outdir}{c}.clean.ATAC_peak.matrix.temp", "w") as o2:

    # 跳过文件头
    next(f1)
    for line in f1:
        line = line.rstrip()
        a = line.split(',')
        chr, start, end = a[0].replace('_', '\t').replace(':', '\t').replace('-', '\t').split('\t', 2)
        
        n = sum(1 for val in a[1:] if float(val) > 0)
        total = len(a) - 1
        per = n / total
        
        if per > 0.1:
            o1.write(f"{chr}\t{start}\t{end}\t{per}\n")
            o2.write(f"{chr}\t{start}\t{end}")
            for val in a[1:]:
                o2.write(f"\t{val}")
            o2.write("\n")

# 使用 sort-bed 命令对临时文件进行排序，并删除临时文件
os.system(f"sort-bed {outdir}{c}.clean.ATAC_peak.temp > {outdir}{c}.clean.ATAC_peak.bed")
os.system(f"sort-bed {outdir}{c}.clean.ATAC_peak.matrix.temp > {outdir}{c}.clean.ATAC_peak.matrix.txt")
os.remove(f"{outdir}{c}.clean.ATAC_peak.temp")
os.remove(f"{outdir}{c}.clean.ATAC_peak.matrix.temp")

# 调用外部的脚本
os.system(f"Rscript 2.R {outdir}{c}.clean.ATAC_peak.bed {reference}")

