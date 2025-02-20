# run_analysis.py
'''import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
import pandas as pd
import sys



# 激活 pandas 和 R 数据框之间的转换
pandas2ri.activate()
ro.r('.libPaths(c("/mnt/kongtianci/R/x86_64-pc-linux-gnu-library/4.2"))')
ro.r('library(WGCNA)')
ro.r('library(gplots)')
ro.r('library(ggplot2)')
ro.r('library(reshape)')
ro.r('library(RColorBrewer)')

# 加载 R 脚本中的函数
ro.r('source("analysis.R")')

# 获取 R 中的函数引用
analyze_expression = ro.globalenv['analyze_expression']

# 解析命令行参数
if len(sys.argv) != 4:
    print("Usage: python run_analysis.py <count_file> <pheno_file> <output_dir>")
    sys.exit(1)

count_file = sys.argv[1]
pheno_file = sys.argv[2]
output_dir = sys.argv[3]

# 加载 count 和 pheno 数据
count = pd.read_csv(count_file, sep='\t')
pheno = pd.read_csv(pheno_file, sep='\t')

# 转换为 R 数据框
r_count = pandas2ri.py2rpy(count)
r_pheno = pandas2ri.py2rpy(pheno)

# 调用 R 函数并传递参数
analyze_expression(r_count, r_pheno, output_dir)

print("分析完成！")
'''

# run_analysis.py
import argparse
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
import pandas as pd

# 激活 pandas 和 R 数据框之间的转换
pandas2ri.activate()



# 加载 R 脚本中的函数
ro.r('source("analysis.R")')
analyze_expression = ro.globalenv['analyze_expression']  # 获取 R 函数引用

# 使用 argparse 解析命令行参数
def parse_arguments():
    parser = argparse.ArgumentParser(description="Run WGCNA analysis with given input files.")
    parser.add_argument("--count_file", "-c", type=str, help="Path to the count matrix file.")
    parser.add_argument("--pheno_file", "-p", type=str, help="Path to the phenotype metadata file.")
    parser.add_argument("--trait", "-t", type=str, help="Relevance of selected trait data.")
    parser.add_argument("--output_dir", "-o", type=str, help="Directory to save the results.")
    return parser.parse_args()

def main():
    # 解析命令行参数
    args = parse_arguments()

    # 设置 R 包路径并加载所需的 R 包
    ro.r('.libPaths(c("/mnt/kongtianci/R/x86_64-pc-linux-gnu-library/4.2"))')
    ro.r('library(WGCNA)')
    ro.r('library(gplots)')
    ro.r('library(ggplot2)')
    ro.r('library(reshape)')
    ro.r('library(RColorBrewer)')
    print(f"Trait passed to R: {args.trait}")
    # 加载 count 和 pheno 数据
    count = pd.read_csv(args.count_file, sep='\t')
    pheno = pd.read_csv(args.pheno_file, sep='\t')

    # 转换为 R 数据框
    r_count = pandas2ri.py2rpy(count)
    r_pheno = pandas2ri.py2rpy(pheno)

    # 调用 R 函数并传递参数
    analyze_expression(r_count, r_pheno, args.trait, args.output_dir)

    print("分析完成！")

if __name__ == "__main__":

    main()
