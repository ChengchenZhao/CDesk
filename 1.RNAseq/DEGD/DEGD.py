import argparse
import os

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


def deseq_analysis(count_file:str, treatment_groups:list, control_groups:list, output_dir:str, **kwargs):
    import omicverse as ov
    # 读取count文件，转换为int格式
    data=pd.read_csv(count_file,index_col=0,sep='\t',header=0)
    data = data.astype(int)
    data_filtered = data[treatment_groups+control_groups]
    data_filtered = data_filtered[data_filtered.sum(axis=1)>=1]
    
    # 生成deg对象，并除重复index
    dds=ov.bulk.pyDEG(data_filtered)
    dds.drop_duplicates_index()
    print('... drop_duplicates_index success')
    
    # 使用deseq2做分析
    result = dds.deg_analysis(treatment_groups,control_groups,
            method='DEseq2')
    
    # 过滤低表达基因
    result=result.loc[result['log2(BaseMean)']>1]
        
    # 参考阈值计算，可以不要
    reference_threshold = dds.foldchange_set(fc_threshold=-1,
                   pval_threshold=0.05,
                   logp_max=6)
    
    # 绘制火山图，并存储
    dds.plot_volcano(title='DEG Analysis',figsize=(16,16),
                 plot_genes_num=8,plot_genes_fontsize=12,)
    
    plt.savefig(os.path.join(output_dir, "volcano_plot.pdf"), format="pdf")
    
    # 绘制热图，并存储
    if "gene_names" in kwargs:
        # 指定基因名字模式
        heatmap_plot_data = data_filtered.loc[kwargs.gene_names.split()]
        return
    else:
        # 指定基因数目模式
        gene_nums = int(kwargs.gene_nums) if "gene_nums" in kwargs else 30
        idx = result.sort_values(by="pvalue").index[:gene_nums]
        heatmap_plot_data = data_filtered.loc[idx]
    
    ### log2plus变换
    log2count = np.log2(heatmap_plot_data+1)
    plt.figure(figsize=(16, 12))  # 设置图像大小
    sns.heatmap(log2count, annot=False, cmap="RdYlBu_r", fmt="g", square=True) # 不显示数值、方格绘图
    plt.savefig(os.path.join(output_dir, "heatmap_plot.pdf"), format="pdf")

# 显示热图
plt.show()
        
def main():
    # 创建 ArgumentParser 对象
    parser = argparse.ArgumentParser(description="Script for differential expression gene detection & analysis")
    
    # 添加 --save_dir / -s 参数（用于指定保存目录）
    parser.add_argument(
        '--output_dir', '-o', 
        type=str, 
        required=True,  # 这个参数是必须的
        help='Directory for saving output'
    )
    
    # 添加 --input / -i 参数（用于接受输入文件）
    parser.add_argument(
        '--input', '-i', 
        type=str, 
        required=True,  # 这个参数也是必须的
        help='Input count matrix'
    )
    
    # 添加 --map / -m 参数（用于指定映射文件）
    parser.add_argument(
        '--species', '-s', 
        type=str, 
        required=False,  # 这个参数也是必须的
        help='Species for choosing .gft for mapping gene id to gene name'
    )
    
    # 添加 --mode 参数，用于选择处理模式
    parser.add_argument(
        '--mode', 
        type=str, 
        required=True,  # 必须提供
        choices=['deseq2', 'gfold', 'adjusted_t'],  # 限制只能是这三个值之一
        help='Mode for DGE, only support "deseq2", "gfold", or "adjusted_t"'
    )
    
    # 添加 --meta 参数，用于选择处理模式
    parser.add_argument(
        '--meta', 
        type=str, 
        required=True,  # 必须提供
        help='Meta file for analysis'
    )
    # 解析命令行参数
    args = parser.parse_args()
    
    # 获取 save_dir 和 map 参数的值
    species = args.species
    count_file = args.input
    output_dir = args.output_dir
    mode = args.mode.lower()
    meta = args.meta
    

    # 检查目录是否存在，不存在则创建
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Make folder: {output_dir}")
    else:
        print(f"Folder exists: {output_dir}")

    # 假设 .gtf 文件存在，输出加载信息
    if species:
        print(f"Species is specified as {species}, GTF file will be used for mapping gene IDs to gene names")
    else:
        print(f"Species is not specified, using gene names in count file")
    
    # 读取meta文件，里面存放分组信息
    if not os.path.exists(meta):
        print(f"Meta file not exists: {meta}")
    with open(meta,'r') as f:
        meta_list = [[s.strip() for s in line.strip().split(":")[1].split(",")] for line in f]
        treatment_groups = meta_list[0]
        control_groups = meta_list[1]

    # 调用mode的算法进行运算        
    if mode=='deseq2':
        deseq_analysis(count_file, treatment_groups, control_groups, output_dir,)
    elif mode=="gfold":
        return
    elif mode=="adjusted_t":
        return
    

if __name__ == '__main__':
    main()
