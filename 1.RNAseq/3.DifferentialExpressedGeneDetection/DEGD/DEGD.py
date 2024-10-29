import argparse
import os

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib import gridspec
import omicverse as ov
import subprocess


def deseqMode(args):
        
    # 获取输入
    species = getGTF(args.species)
    count_file = args.input
    output_dir = args.output_dir
    meta_file = args.meta_file
    plot = args.plot
    goi = args.gene_of_interest.split(',')
    top_genes = int(args.top_genes)
    
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
    
    # 把count和fpkm等矩阵存储到输出目录
    count_df = pd.read_csv(count_file, index_col=0, sep='\t', header=0)
    count_df.to_csv(os.path.join(output_dir,'count.csv'),sep='\t')
    # 算fpkm等等...并存储
    #
    
    
    # 读取meta文件，里面存放分组信息，并分析
    metas = readMetaFile(meta_file)
    if plot:
        draw_args = vars(args).copy()
        draw_args['input_dir']=args.output_dir

    for meta in metas:
        result = deseqAnalysis(count_file, species, meta, output_dir, plot=plot)
        gene_list = list(set(result.sort_values(by="qvalue").index[:top_genes]) | set(goi))
        # 如果要画图，那么就画图
        if plot:
            drawHeatmapPlot(gene_list=gene_list, meta=meta, **draw_args)
            drawVocanoPlot(meta=meta, **draw_args)
    return

def deseqAnalysis(count_file:str, species:str, meta:tuple, output_dir:str, plot:bool, **kwargs):
    
    # ---解析meta文件-----------------------------------------
    (name, treatment_groups, control_groups) = meta
    
    # ---读取count文件，转换为int格式--------------------------
    data=pd.read_csv(count_file,index_col=0,sep='\t',header=0)
    data = data.astype(int)
    data_filtered = data[treatment_groups+control_groups]
    data_filtered = data_filtered[data_filtered.sum(axis=1)>=1]
    
    # ---deseq2做分析-----------------------------------------
    
    # 生成deg对象，并除重复index
    dds=ov.bulk.pyDEG(data_filtered)
    dds.drop_duplicates_index()
    print('... drop_duplicates_index success')
    
    # 调用desseq2
    result = dds.deg_analysis(treatment_groups,control_groups,
            method='DEseq2')
    
    # 过滤低表达基因
    result = result.loc[result['log2(BaseMean)']>1]
    
    # 结果存储到本地
    result.to_csv(os.path.join(output_dir,name+'.csv'),sep='\t')
        
    # 参考阈值计算
    # reference_threshold = dds.foldchange_set(fc_threshold=-1,
    #                pval_threshold=0.05,
    #                logp_max=6)
    
    
    # ---绘制火山图，并存储-------------------------------------
    # dds.plot_volcano(title='DEG Analysis',figsize=(16,16),
    #              plot_genes_num=8,plot_genes_fontsize=12,)
    
    # plt.savefig(os.path.join(output_dir, "volcano_plot.pdf"), format="pdf")
    
    return result

def TMode(args):
    
    # 获取输入
    species = getGTF(args.species)
    count_file = args.input
    output_dir = args.output_dir
    meta_file = args.meta_file
    plot = args.plot
    goi = args.gene_of_interest.split(',')
    top_genes = int(args.top_genes)
    
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
    
    # 把count和fpkm等矩阵存储到输出目录
    count_df = pd.read_csv(count_file, index_col=0, sep='\t', header=0)
    count_df.to_csv(os.path.join(output_dir,'count.csv'),sep='\t')
    # 算fpkm等等...并存储
    #
    
    
    # 读取meta文件，里面存放分组信息，并分析
    metas = readMetaFile(meta_file)
    if plot:
        draw_args = vars(args).copy()
        draw_args['input_dir']=args.output_dir

    for meta in metas:
        result = TAnalysis(count_file, species, meta, output_dir, plot=plot)
        gene_list = list(set(result.sort_values(by="qvalue").index[:top_genes]) | set(goi))
        # 如果要画图，那么就画图
        if plot:
            drawHeatmapPlot(gene_list=gene_list, meta=meta, **draw_args)
            drawVocanoPlot(meta=meta, **draw_args)
    return

def TAnalysis(count_file:str, species:str, meta:tuple, output_dir:str, plot:bool, **kwargs):
    
    # ---解析meta文件-----------------------------------------
    (name, treatment_groups, control_groups) = meta
    
    # ---读取count文件，转换为int格式--------------------------
    data=pd.read_csv(count_file,index_col=0,sep='\t',header=0)
    data = data.astype(int)
    data_filtered = data[treatment_groups+control_groups]
    data_filtered = data_filtered[data_filtered.sum(axis=1)>=1]
    
    # ---deseq2做分析-----------------------------------------
    
    # 生成deg对象，并除重复index
    dds=ov.bulk.pyDEG(data_filtered)
    dds.drop_duplicates_index()
    print('... drop_duplicates_index success')
    
    # 调用desseq2
    result = dds.deg_analysis(treatment_groups,control_groups,
            method='ttest',multipletests_method='bonferroni')
    
    # 过滤低表达基因
    result = result.loc[result['log2(BaseMean)']>1]
    
    # 结果存储到本地
    result.to_csv(os.path.join(output_dir,name+'.csv'),sep='\t')
        
    return result

def gfoldMode(args):
    
    # 获取参数
    
    species = getGTF(args.species)
    input_dir = args.input_dir
    output_dir = args.output_dir
    meta_file = args.meta_file
    plot = args.plot
    goi = args.gene_of_interest.split(',')
    top_genes = int(args.top_genes)
    
    # 将输入目录中的所有bam文件存入list
    bam_files = [os.path.join(input_dir,f) for f in os.listdir(input_dir) if f.endswith('.bam')]

    # 检查输出目录是否存在，不存在则创建
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Make folder: {output_dir}")
    else:
        print(f"Folder already exists: {output_dir}")
    if not os.path.exists(os.path.join(output_dir,"readCnt")):    
        os.makedirs(os.path.join(output_dir,"readCnt"))
    if not os.path.exists(os.path.join(output_dir,"diff")):
        os.makedirs(os.path.join(output_dir,"diff"))

    # 假设 .gtf 文件存在，输出加载信息
    if species:
        print(f"Species is specified as {species}, GTF file will be used for mapping gene IDs to gene names")
    else:
        print(f"Species file not exist!")
        return
    
    # 读取meta文件，里面存放分组信息   
    metas = readMetaFile(meta_file)
    
    # bam文件转count matrix
    gfoldBam2Cnt(bam_files=bam_files, species=species, output_dir=output_dir)
    
    if plot:
        draw_args = vars(args).copy()
        draw_args['input_dir']=args.output_dir
    
    # 做差异分析
    for meta in metas:
        result = gfoldAnalysis(meta=meta, **vars(args))
        gene_list = list(set(result['GFOLD'].abs().nlargest(top_genes).index) | set(goi))
        
    # 画图
        if plot:
            drawHeatmapPlot(gene_list=gene_list, meta=meta, **draw_args)
            drawVocanoPlot(meta=meta,fc_type='GFOLD', **draw_args)    
    
    return

def gfoldBam2Cnt(bam_files:list, species:str, output_dir:str):
    # GFOLD方法需要大量调用命令行中的GFOLD
    
    # ---首先将bam转换为 read_cnt文件----------------------------------------

    processes = []
    for file in bam_files:
        name = file.split("/")[-1].split(".")[0]
        # command = f"samtools view -@ 20 {output_path}/bam/{file} | gfold count -ann ${gtf} -tag stdin -o ${name}.read_cnt"
        command = f"samtools view -@ 20 {file} | gfold count -ann {species} -tag stdin -o {os.path.join(output_dir, 'readCnt', name + '.read_cnt')}"

        processes.append(subprocess.Popen(command, shell=True))
    
    # ---等待所有进程完成，将所有结果整合入count_matrix方便后续分析--------------
    
    for p in processes:
        p.wait()
    
    # 初始化一个空的 DataFrame
    combined_count_df = pd.DataFrame()
    combined_fpkm_df = pd.DataFrame()
    for file in os.listdir(os.path.join(output_dir, 'readCnt')):
        if file.endswith('.read_cnt'):
            name = file[:-9]
            
            # 读取文件
            df = pd.read_csv(os.path.join(output_dir, 'readCnt',file), header=None, index_col=0, sep='\t')
            count_df = df[2]
            count_df.name = name
            fpkm_df = df[4]
            fpkm_df.name = name
            
            # 合并到总的 DataFrame
            combined_count_df = pd.concat([combined_count_df, count_df], axis=1)
            combined_fpkm_df = pd.concat([combined_fpkm_df, fpkm_df], axis=1)
    
    # 存储到output
    combined_count_df.to_csv(os.path.join(output_dir,'count.csv'),sep='\t')
    combined_fpkm_df.to_csv(os.path.join(output_dir,'fpkm.csv'),sep='\t')

def gfoldAnalysis(meta:tuple, species:str, output_dir:str, **kwargs):
    
    (name, treatment_groups, control_groups) = meta
    # ---基于meta分组信息进行差异分析------------------------------------------
    os.chdir(os.path.join(output_dir,"readCnt"))
    command = "gfold diff \
    -s1 {} \
    -s2 {} \
    -suf .read_cnt \
    -o ../diff/{}.diff \
    > ../diff/{}.gfold.log".format(
        ",".join(treatment_groups),
        ",".join(control_groups),
        name, name
        )
    subprocess.run(command, shell=True)
    
    # ---统计结果生成result文档，方便后续调用绘图------------------------------
    result = deseqAnalysis(os.path.join(output_dir,'count.csv'), species, meta, output_dir,plot=False)

    df = pd.read_csv("../diff/{}.diff".format(name), sep="\t", comment="#",header=None, index_col=1)
    df = df[[2,3,4]]
    df.columns = ["GFOLD", "E-FDR", "log2FoldChange"]
    result = pd.concat([result,df], join='inner', axis=1)
    result.to_csv(os.path.join(output_dir,name+'.csv'),sep='\t')
    
    # ---统计GFOLD值特异的基因------------------------------------------------
    # os.chdir("../diff")

    # subprocess.run('grep -v "#" TreatVSControl.diff | awk -v "OFS=\\t" \'{if ($3 > 1.5 || $3 < -1.5) print $1}\' > TreatVSControl_cutoff_15.deg', shell=True)
    # subprocess.run('grep -v "#" TreatVSControl.diff | awk -v "OFS=\\t" \'{if ($3 > 1 || $3 < -1) print $1}\' > TreatVSControl_cutoff_1.deg', shell=True)
    # subprocess.run('grep -v "#" TreatVSControl.diff | awk -v "OFS=\\t" \'{if ($3 > 0.5 || $3 < -0.5) print $1}\' > TreatVSControl_cutoff_05.deg', shell=True)
    
    return result

def drawMode(args):
    
    analysis_path = args.input_dir
    method = args.analysis_method
    output_dir = args.output_dir
    heatmap_types = args.heatmap_types
    volcano_types = args.volcano_types
    meta_file = args.meta_file
    
    gene_nums = args.gene_nums
    gene_names = args.gene_names
    
    if os.path.exists(analysis_path):
        raise ValueError('Input path not exist! Please check the path')


    if method=='deseq2':
        # 生成要做绘图展示的gene list
        
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
        gene_list=[]
        
        # 针对所有输入的参数绘图
        drawHeatmapPlot(analysis_path=analysis_path, 
                        heatmap_types=heatmap_types, 
                        gene_list=gene_list, 
                        meta_file=meta_file, 
                        output_dir=output_dir)
        
    elif method=='gfold':
        return
    pass

def drawHeatmapPlot(input_dir:str, heatmap_types:str, gene_list:list, meta:tuple, output_dir:str,**kwargs):
    
    heatmap_types = heatmap_types.split(',')
    
    # 动态设置子图布局
    num_plots = len(heatmap_types)
    cols = 3  # 设置列数
    rows = (num_plots + cols - 1) // cols  # 计算行数

    fig = plt.figure(figsize=(8 * cols, 16 * rows))
    gs = gridspec.GridSpec(rows, cols)
    
    (name, treatment_groups, control_groups) = meta
    # 绘制每个 heatmap
    for i, heatmap_type in enumerate(heatmap_types):
        
        data_path = os.path.join(input_dir,heatmap_type+'.csv')
        if not os.path.exists(data_path):
            raise ValueError(f'{heatmap_type} data not exists!')
        data = pd.read_csv(data_path, sep="\t", header=0, index_col=0)
        data = data[treatment_groups+control_groups].loc[gene_list]
        
        ### log2plus变换
        data = np.log2(data+1) 

        ax = fig.add_subplot(gs[i // cols, i % cols])  # 计算子图的位置
        sns.heatmap(data, ax=ax, annot=False, cmap="RdYlBu_r", fmt="g", square=True)
        ax.set_title(f'{heatmap_type} Heatmap')

    # 去掉多余的子图
    #for j in range(num_plots, rows * cols):
    #    if j < len(fig.axes):  # 检查索引是否在范围内
    #        fig.delaxes(fig.axes[j])

    # 显示图形
    plt.tight_layout()

    plt.savefig(os.path.join(output_dir, name+"_heatmap.pdf"), format="pdf")

def drawVocanoPlot(input_dir:str, meta:tuple, output_dir:str, vocano_types:str='qvalue', fc_type='log2FC',
                   pval_threshold:float=0.01, fc_max:float=1.5, fc_min:float=-1.5, fc_threshold=-1,
                   logp_max:int=6,**kwargs):
    
    vocano_types = vocano_types.split(',')
    
        # 动态设置子图布局
    num_plots = len(vocano_types)
    cols = 2  # 设置列数
    rows = (num_plots + cols - 1) // cols  # 计算行数

    fig = plt.figure(figsize=(16 * cols, 16 * rows))
    gs = gridspec.GridSpec(rows, cols)
    
    (name, treatment_groups, control_groups) = meta
    # 绘制每个 heatmap
    fold_threshold = 0
    for i, vocano_type in enumerate(vocano_types):
        result = pd.read_csv(os.path.join(input_dir,name+'.csv'),header=0,index_col=0,sep='\t')
        if fc_threshold==-1:
            foldp=np.histogram(result[fc_type].dropna())
            foldchange=(foldp[1][np.where(foldp[1]>0)[0][fold_threshold]]+foldp[1][np.where(foldp[1]>0)[0][fold_threshold+1]])/2
        else:
            foldchange=float(fc_threshold)
        print('... Fold change threshold: %s'%foldchange)
        fc_max,fc_min=foldchange,0-foldchange
        pval_threshold=pval_threshold
        result['sig']='normal'
    
        result.loc[((result[fc_type]>fc_max)&(result[vocano_type]<pval_threshold)),'sig']='up'
        result.loc[((result[fc_type]<fc_min)&(result[vocano_type]<pval_threshold)),'sig']='down'
        result.loc[result['-log(qvalue)']>logp_max,'-log(qvalue)']=logp_max
        
        ax = fig.add_subplot(gs[i // cols, i % cols])  # 计算子图的位置
        
        ov.pl.volcano(result,ax=ax, pval_name='qvalue',fc_name=fc_type,
                     pval_threshold=pval_threshold,fc_max=fc_max,fc_min=fc_min,
                      pval_max=6,FC_max=10,
                    figsize=(8,8),title='DEGs in Bulk',titlefont={'weight':'normal','size':14,},
                     up_color='#e25d5d',down_color='#7388c1',normal_color='#d7d7d7',
                     up_fontcolor='#e25d5d',down_fontcolor='#7388c1',normal_fontcolor='#d7d7d7',
                     legend_bbox=(0.8, -0.2),legend_ncol=2,legend_fontsize=12,
                     plot_genes=None,plot_genes_num=10,plot_genes_fontsize=11,
                     ticks_fontsize=12,)

        ax.set_title(f'{vocano_type} vocano')
        ax.set_ylabel(vocano_type)                                    
        ax.set_xlabel(fc_type) 

    # 去掉多余的子图
    #for j in range(num_plots, rows * cols):
    #    if j < len(fig.axes):  # 检查索引是否在范围内
    #        fig.delaxes(fig.axes[j])

    # 显示图形
    plt.tight_layout()

    plt.savefig(os.path.join(output_dir, name+"_vocano.pdf"), format="pdf")
 
def readMetaFile(meta):
    if not os.path.exists(meta):
        raise ValueError(f"Meta file not exists: {meta}")
    df = pd.read_excel(meta,header=0,index_col=0)
    result = []
    for column in df.columns[1:]:
        # 获取该列的值
        values = df[column]
        
        # 生成实验组和对照组的样本名列表
        treatment_groups = values.index[values==1].to_list()
        control_groups = values.index[values==-1].to_list()
        
        # 创建元组并添加到结果列表
        result.append((str(column), treatment_groups, control_groups))
    return result

def getGTF(species):
    if species=="mouse":
        gtf = os.path.join(os.path.dirname(os.path.abspath(__file__)),os.pardir,"gtf","mm10.ncbiRefSeq.WithUCSC.gtf")
    return gtf

def addCommonArguments(parser, func):
    
    # 添加 --save_dir / -s 参数（用于指定保存目录）
    parser.add_argument(
        '--output_dir', '-o', 
        type=str, 
        required=True,  # 这个参数是必须的
        help='Directory for saving output'
    )
    
    # 添加 --species / -s 参数（用于指定映射文件）
    parser.add_argument(
        '--species', '-s', 
        type=str, 
        default="mouse",  # 默认为小鼠
        help='Species for choosing .gtf for mapping gene id to gene name'
    )
    
    # 添加 --meta 参数，用于提供分组信息
    parser.add_argument(
        '--meta_file', '-m',
        type=str, 
        required=True,  # 必须提供
        help='Group informations for analysis'
    )
    
    # 添加 --plot 参数，决定是否要绘图
    parser.add_argument(
        '--plot', '-p',
        action='store_true', 
        help='Whether draw plots or not'
    )
    
    # 添加 --heatmap_types 参数，决定是否要绘图
    parser.add_argument(
        '--heatmap_types', '-H',
        default='count', 
        help='Data types for drawing heatmap plots'
    )
    
    # 添加 --vocano_types 参数，决定是否要绘图
    parser.add_argument(
        '--vocano_types', '-v',
        default='qvalue', 
        help='Data types for drawing vocano plots'
    )
    
    # 添加 --gene_of_interest 参数，决定是否要绘图
    parser.add_argument(
        '--gene_of_interest', '--goi', '-g',
        default='', 
        help='Interest gene for plot heatmap'
    )
    # 添加 --top_genes 参数，决定是否要绘图
    parser.add_argument(
        '--top_genes', '-t',
        default=30, 
        help='Number of top differential expression genes for plot heatmap'
    )
    
    # 添加 --fc_thresold 参数，决定是否要绘图
    parser.add_argument(
        '--fc_threshold',
        default=1.5, 
        help='fc_threshold'
    )
    
    # 指定默认分析函数
    parser.set_defaults(func=func)
      
def main():
    # 创建 ArgumentParser 对象
    parser = argparse.ArgumentParser(description="Script for differential expression gene detection & analysis")
    
    # 创建子解析器
    subparsers = parser.add_subparsers(dest='mode', required=True, help="Choose a mode to run the script")


    # -------------------------------------------
    # Mode DESeq2 子解析器
    # -------------------------------------------
    
    parser_deseq2 = subparsers.add_parser('deseq2', help="Run by deseq2")
    
    # 添加 --input / -i 参数（用于接受输入文件）
    parser_deseq2.add_argument(
        '--input', '-i', 
        type=str, 
        required=True,  # 这个参数也是必须的
        help='Input count matrix for DESeq2'
    )

    addCommonArguments(parser_deseq2, func=deseqMode)


    # -------------------------------------------
    # Mode GFOLD 子解析器
    # -------------------------------------------
    parser_gfold = subparsers.add_parser('gfold', help="Run by gfold")
    
    parser_gfold.add_argument(
        '--input_dir', '-i', 
        type=str,
        required=True, 
        help="Input bam file directory for GFOLD")
    
    addCommonArguments(parser_gfold,func=gfoldMode)
    
    
    # -------------------------------------------
    # Mode adjusted_t 子解析器
    # -------------------------------------------
    
    parser_adjusted_t = subparsers.add_parser('adjusted_t', help="Run by adjusted_t")
    
    # 添加 --input / -i 参数（用于接受输入文件）
    parser_adjusted_t.add_argument(
        '--input', '-i', 
        type=str, 
        required=True,  # 这个参数也是必须的
        help='Input count matrix for adjusted_t'
    )

    addCommonArguments(parser_adjusted_t, func=TMode)
    
    
    # -------------------------------------------
    # Mode plot 子解析器
    # -------------------------------------------
    
    parser_plot = subparsers.add_parser('plot', help="Draw plots from upstream output")
    
    # 添加 --input_dir
    parser_plot.add_argument(
        '--input_dir', '-i', 
        type=str, 
        required=True,  # 这个参数是必须的
        help='Directory for output of upstream differential analysis'
    )
    
    # 添加 --output_dir
    parser_plot.add_argument(
        '--output_dir', '-o', 
        type=str, 
        required=True,  # 这个参数是必须的
        help='Directory for saving output'
    )
    
    # 添加 --analysis_method
    parser_plot.add_argument(
        '--analysis_method', '-a', 
        type=str, 
        required=True,  # 这个参数是必须的
        help='Method used for differential analysis'
    )
    
    parser_plot.set_defaults(func=drawMode)
    
    
    # 解析命令行参数
    args = parser.parse_args()

    # 检查是否指定了子命令
    if args.mode is None:
        parser.print_help()  # 打印帮助信息
        print("\nError: Must specify a mode ('deseq2' , 'gfold' or 'adjusted_t').")
    
    
    # ---调用子func-------------------------------------
    args.func(args)   

if __name__ == '__main__':
    main()
