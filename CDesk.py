import os
import argparse
import sys
import subprocess

script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)

# bulkRNA 1. preprocessing function
def bulkRNA_preprocess(args):
    bash_script = './1.RNAseq/1.PreProcessPipeline/BulkRNAseqData.sh'
    cmd = [bash_script, '-s', args.s, '-i', args.i, '-o', args.o, '-p', str(args.t), '-l', str(args.l)]
    subprocess.run(cmd,check=True)

    if args.bw:
        bash_script = './1.RNAseq/1.PreProcessPipeline/BulkRNAseqBW.sh'
        cmd = [bash_script, '-i', f"{args.o}/Bam", '-o', f"{args.o}/BW", '-p', str(args.t)]
        subprocess.run(cmd, check=True)

    if args.te:
        bash_script = './1.RNAseq/1.PreProcessPipeline/BulkRNAseqTE.sh'
        cmd = [bash_script, '-i', f"{args.o}/Bam", '-o', f"{args.o}/TE", '-p', str(args.t)]
        subprocess.run(cmd, check=True)

# bulkRNA 2. QC
def bulkRNA_QC(args):
    bash_script = './1.RNAseq/2.QuanlityControl/BulkRNAseqRep.sh'
    cmd = [ 
    bash_script, 
    '-i', args.i,
    '-r', args.r,
    '-o', args.o]
    subprocess.run(cmd, check=True)

# bulkRNA 3. DEG
def bulkRNA_DEG(args):
    if args.DEG_mode == 'deseq2':
        python_script = './1.RNAseq/3.DifferentialExpressedGeneDetection/DEGD.py'
        if args.p:
            cmd = [
                "python", python_script, "deseq2",
                "--input",args.i, "--output_dir",args.o,
                "--meta",args.m,"-p",
                "--heatmap_types",args.H,"--volcano_types",args.v,
                "--gene_of_interest",args.g,"--top_genes",str(args.t),
                "--fc_threshold",str(args.fc)
            ]
            subprocess.run(cmd, check=True)
        else:
            cmd = [
                "python", python_script, "deseq2",
                "--input",args.i, "--output_dir",args.o,
                "--meta",args.m
            ]
            subprocess.run(cmd, check=True)
    elif args.DEG_mode == 'adjusted_t':
        python_script = './1.RNAseq/3.DifferentialExpressedGeneDetection/DEGD.py'
        if args.p:
            cmd = [
                "python", python_script, "adjusted_t",
                "--input",args.i, "--output_dir",args.o,
                "--meta",args.m,"-p",
                "--heatmap_types",args.H,"--volcano_types",args.v,
                "--gene_of_interest",args.g,"--top_genes",str(args.t),
                "--fc_threshold",str(args.fc)
            ]
            subprocess.run(cmd, check=True)
        else:
            cmd = [
                "python", python_script, "adjusted_t",
                "--input",args.i, "--output_dir",args.o,
                "--meta",args.m
            ]
            subprocess.run(cmd, check=True)
    elif args.DEG_mode == 'gfold':
        python_script = './1.RNAseq/3.DifferentialExpressedGeneDetection/DEGD.py'
        if args.p:
            cmd = [
                "python", python_script, "gfold",
                "--input",args.i, "--output_dir",args.o,
                "--meta",args.m,"-p","--species",args.s,
                "--heatmap_types",args.H,"--volcano_types",args.v,
                "--gene_of_interest",args.g,"--top_genes",str(args.t),
                "--fc_threshold",str(args.fc)
            ]
            subprocess.run(cmd, check=True)
        else:
            cmd = [
                "python", python_script, "gfold",
                "--input",args.i, "--output_dir",args.o,
                "--meta",args.m,"--species",args.s
            ]
            subprocess.run(cmd, check=True)

# bulkRNA 5. GO analysis
def bulkRNA_GO(args):
    if args.GO_mode == 'custom':
        R_script = './1.RNAseq/5.GeneOntologyAnalysis/enrichment_analysis_customer.R'
        cmd = [
            "Rscript", R_script, 
            args.i, args.c, args.o
        ]
        subprocess.run(cmd, check=True)
    elif args.GO_mode == 'GO':
        R_script = './1.RNAseq/5.GeneOntologyAnalysis/enrichment_analysis.R'
        cmd = [
            "Rscript", R_script, 
            args.i, 'GO', args.s, str(args.t), args.o
        ]
        subprocess.run(cmd, check=True)
    elif args.GO_mode == 'KEGG':
        R_script = './1.RNAseq/5.GeneOntologyAnalysis/enrichment_analysis.R'
        cmd = [
            "Rscript", R_script, 
            args.i, 'KEGG', args.s, args.o
        ]
        subprocess.run(cmd, check=True)
    elif args.GO_mode == 'plot':
        R_script = './1.RNAseq/5.GeneOntologyAnalysis/plot_enrichment.R'
        cmd = [
            "Rscript", R_script, 
            args.i, args.color_high, args.color_low, args.o
        ]
        subprocess.run(cmd, check=True)

# bulkRNA 6. Similarity Analysis
def bulkRNA_Similarity(args):
    bash_script = './1.RNAseq/6.SimilarityAnalysis/BulkRNAseqSim.sh'
    cmd = [bash_script, '-s', args.s, '-i', args.i, '-o', args.o, '-p', args.t, '-l', args.l]
    subprocess.run(cmd, check=True)

# bulkRNA 7. WGCNA
def bulkRNA_WGCNA(args):
    python_script = './1.RNAseq/7.GeneClustering/run_analysis.py'
    cmd = ['python', python_script, '--c', args.count_file, '--p', args.pheno_file, '--t', args.trait, '--o', args.o]
    subprocess.run(cmd, check=True)

# bulkRNA 8. AlternativeSplicingDetection
def bulkRNA_AlternativeSpliceDetection(args):
    if args.Splice_mode == 'detect':
        bash_script = './1.RNAseq/8.AlternativeSplicingDetection/alternative_splicing.sh'
        b1_file = "b1.txt"
        b2_file = "b2.txt"
        with open(b1_file, "w") as f1:
            f1.write(args.b1)
        with open(b2_file, "w") as f2:
            f2.write(args.b2)
        cmd = [bash_script,'--b1',b1_file,'--b2',b2_file,'--gtf',args.gtf,'--od',args.o]
        subprocess.run(cmd,check=True)
        os.remove(b1_file)
        os.remove(b2_file)
    elif args.Splice_mode == 'draw':
        bash_script = './1.RNAseq/8.AlternativeSplicingDetection/draw_sashimiplot.sh'
        #cmd = (f"{bash_script} --b1 {args.b1} --b2 {args.b2} --gff {args.gff} --r {args.interval} --g {args.group} --od {args.o}")
        #os.system(cmd)
        cmd = [bash_script,'--b1',args.b1,'--b2',args.b2,'--gff',args.gff,'--r',args.interval,'--g',args.group,'--od',args.o]
        subprocess.run(cmd,check=True)

# scRNA 1. preprocess
def scRNA_preprocess(args):
    python_script = './2.scRNA/1.preprocess/preprocess.py'
    if args.scRNA_preprocess_mode == 'cellranger':
        cmd = ["python", python_script,"cellranger","-i",args.i,"-o",args.o,"-s",args.s,'-t',str(args.t)]
        subprocess.run(cmd,check=True)
    elif args.scRNA_preprocess_mode == 'drseq':
        cmd = ["python", python_script,"drseq","-i",args.i,"-o",args.o,"-s",args.s,'-t',str(args.t)]
        subprocess.run(cmd,check=True)
    elif args.scRNA_preprocess_mode == 'smartseq':
        bash_script = './1.RNAseq/1.PreProcessPipeline/BulkRNAseqData.sh'
        cmd = [bash_script, '-s', args.s, '-i', args.i, '-o', args.o, '-p', str(args.t), '-l', str(args.l)]
        subprocess.run(cmd,check=True)
        if args.bw:
            bash_script = './1.RNAseq/1.PreProcessPipeline/BulkRNAseqBW.sh'
            cmd = [bash_script, '-i', f"{args.o}/Bam", '-o', f"{args.o}/BW", '-p', str(args.t)]
            subprocess.run(cmd, check=True)
        if args.te:
            bash_script = './1.RNAseq/1.PreProcessPipeline/BulkRNAseqTE.sh'
            cmd = [bash_script, '-i', f"{args.o}/Bam", '-o', f"{args.o}/TE", '-p', str(args.t)]
            subprocess.run(cmd, check=True)
    elif args.scRNA_preprocess_mode == 'singleron':
        cmd = ["python", python_script,"singleron","-i",args.i,"-o",args.o,"-s",args.s,"-t",str(args.t)]
        subprocess.run(cmd, check=True)

# scRNA 3. cluster
def scRNA_cluster(args):
    python_script = './2.scRNA/3.clustering/clustering.py'
    R_script = './2.scRNA/3.clustering/clustering.R'
    cmd = [
    'python', python_script,
    '-i', args.i,  
    '-o', args.o,  
    '-n', args.n,  
    '--min_cells', str(args.min_cells),  
    '--min_features', str(args.min_features),  
    '--nFeature_RNA_min', str(args.nFeature_RNA_min),  
    '--nFeature_RNA_max', str(args.nFeature_RNA_max),  
    '--mt_percent', str(args.mt_percent), 
    '--variable_features', str(args.variable_features),  
    '--dim_prefer', str(args.dim_prefer), 
    '--res', str(args.res) 
    ]
    subprocess.run(cmd, check=True)
    cmd = [
    'Rscript', R_script,
    args.i, 
    args.o,  
    args.n,  
    str(args.min_cells),  
    str(args.min_features),
    str(args.nFeature_RNA_min), 
    str(args.nFeature_RNA_max),  
    str(args.mt_percent),  
    str(args.variable_features),  
    str(args.dim_prefer),  
    str(args.res)  
    ]
    subprocess.run(cmd, check=True)

# scRNA: 7. trajectory
def scRNA_trajectory(args):
    if args.scRNA_trajectory_mode == 'monocle':
        R_script = './2.scRNA/7.trajectory/monocle3.R'
        cmd = ['conda','run','-n','monocle','Rscript',R_script,args.i,args.o,args.m,args.s]
        subprocess.run(cmd, check=True)
    elif args.scRNA_trajectory_mode == 'diffusion':
        R_script = './2.scRNA/7.trajectory/diffusionmap.R'
        cmd = ['conda','run','-n','monocle','Rscript',R_script,args.i,args.o,args.m,args.s]
        subprocess.run(cmd, check=True)
    elif args.scRNA_trajectory_mode == 'velocity':
        python_script = './2.scRNA/7.trajectory/velocity.py'
        cmd = ['conda','run','-n','velocyto','python',python_script,'--cellranger_output',args.cellranger,
               '--rds_file',args.rds,'--gtf',args.gtf,'--genes',args.genes,'-o',args.o]
        subprocess.run(cmd, check=True)
    elif args.scRNA_trajectory_mode == 'cstreet':
        python_script = './2.scRNA/7.trajectory/cstreet.py'
        cmd = ['conda','run','-n','cstreet','python',python_script,'-i',args.i,'-s',args.s,'-n',args.n,'-o',args.o]
        subprocess.run(cmd, check=True)

# tools 1. bw
def bw(args):
    bash_script = './1.RNAseq/1.PreProcessPipeline/BulkRNAseqBW.sh'
    cmd = (f"{bash_script} -i {args.i} -o {args.o} -p {args.p}")
    os.system(cmd)

# tools 2. te
def te(args):
    bash_script = './1.RNAseq/1.PreProcessPipeline/BulkRNAseqTE.sh'
    cmd = (f"{bash_script} -i {args.i} -o {args.o} -p {args.p}")
    os.system(cmd)

def main():
    # CDesk ArgumentParser Object
    parser = argparse.ArgumentParser(prog='CDesk', description="CDesk multiomics pipeline")
    # 创建二级子命令
    subparsers = parser.add_subparsers(dest='omics', required=True, help="Choose an omic to run the script")

    # bulkRNA
    # 创建 bulkRNA 子命令的子解析器
    bulkRNA_parser = subparsers.add_parser('bulkRNA', help='bulkRNA-seq pipeline')
    bulkRNA_subparsers = bulkRNA_parser.add_subparsers(help='bulkRNA-seq pipeline', title='bulkRNA-seq pipeline', dest='bulkRNA_command')

    # bulkRNA: 1. preprocess
    # 添加 bulkRNA preprocessing 子命令
    bulkPreprocessing_parser = bulkRNA_subparsers.add_parser('preprocess', help='bulkRNA-seq preprocessing')
    # Required parameters
    bulkPreprocessing_parser.add_argument('-i', required=True, help='Specify the absolute directory of input file (fastq.gz), the ending does not require a /')
    bulkPreprocessing_parser.add_argument('-o', required=True, help='Specify the absolute directory of output file, the ending does not require a /')
    bulkPreprocessing_parser.add_argument('-s', required=True, help='Specify the species',
                                      choices=["mm9", "mm10","pig",'human'])
    # Optional parameters
    bulkPreprocessing_parser.add_argument('-bw', action='store_true', help='If specified, perform a bam -> bigwig file transfer')
    bulkPreprocessing_parser.add_argument('-te', action='store_true', help='If specified, perform a TE expression analysis')
    bulkPreprocessing_parser.add_argument('-t', help='Specify number of threads (default is 8)',default=8,type=int)
    bulkPreprocessing_parser.add_argument('-l', help='Single sequencing, 2:Pair sequencing (default is 2)',default=2)

    # bulkRNA: 2. QC
    # 添加 bulkRNA QC 子命令
    bulkQC_parser = bulkRNA_subparsers.add_parser('QC', help='bulkRNA-seq QC')
    # Required parameters
    bulkQC_parser.add_argument('-i', required=True, help='Specify the input file (.csv), column name as sample, row name as gene')
    bulkQC_parser.add_argument('-o', required=True, help='Specify the absolute directory of output file, the ending does not require a /')
    # Opitional parameters
    bulkQC_parser.add_argument('-r', help='Specify the symbol representing repeated experiments.',default='Rep')

    # bulkRNA: 3. DEG
    bulkDEGD_parser = bulkRNA_subparsers.add_parser('DEG', help='bulkRNA-seq differential gene expression analysis')
    bulkDEGD_subparser = bulkDEGD_parser.add_subparsers(dest='DEG_mode', required=True, help="Choose a mode to run the script")
    # DEseq2
    bulkDEGD_deseq2 = bulkDEGD_subparser.add_parser('deseq2',help='Perform DEseq2 differential gene expression analysis')
    bulkDEGD_deseq2.add_argument('-i',required=True,help='Input count matrix for DESeq2')
    bulkDEGD_deseq2.add_argument('-o',required=True,help='Directory for saving output')
    bulkDEGD_deseq2.add_argument('-m',required=True,help='Group informations for analysis')
    bulkDEGD_deseq2.add_argument('-p',help='Whether to draw plot or not',action='store_true')
    bulkDEGD_deseq2.add_argument('-H',help='Data types for drawing heatmap plots',default='count')    
    bulkDEGD_deseq2.add_argument('-v',help='Data types for drawing volcano plots',default='qvalue')
    bulkDEGD_deseq2.add_argument('-g',help='Interest gene for plot heatmap,split by ","',default='')
    bulkDEGD_deseq2.add_argument('-t',help='Number of top differential expression genes for plot heatmap',default=30,type=int)
    bulkDEGD_deseq2.add_argument('-fc',help='fc_threshold for vocano plot',default=1.5)
    # adjusted_t
    bulkDEGD_adjust_t = bulkDEGD_subparser.add_parser('adjusted_t',help='Perform adjusted_t differential gene expression analysis')
    bulkDEGD_adjust_t.add_argument('-i',required=True,help='Input count matrix for DESeq2')
    bulkDEGD_adjust_t.add_argument('-o',required=True,help='Directory for saving output')
    bulkDEGD_adjust_t.add_argument('-m',required=True,help='Group informations for analysis')
    bulkDEGD_adjust_t.add_argument('-p',help='Whether to draw plot or not',action='store_true')
    bulkDEGD_adjust_t.add_argument('-H',help='Data types for drawing heatmap plots',default='count')
    bulkDEGD_adjust_t.add_argument('-v',help='Data types for drawing volcano plots',default='qvalue')
    bulkDEGD_adjust_t.add_argument('-g',help='Interest gene for plot heatmap,split by ","',default='')
    bulkDEGD_adjust_t.add_argument('-t',help='Number of top differential expression genes for plot heatmap',default=30,type=int)
    bulkDEGD_adjust_t.add_argument('-fc',help='fc_threshold for vocano plot',default=1.5)
    # gfold
    bulkDEGD_adjust_t = bulkDEGD_subparser.add_parser('gfold',help='Perform adjusted_t differential gene expression analysis')
    bulkDEGD_adjust_t.add_argument('-i',required=True,help='Input count matrix for DESeq2')
    bulkDEGD_adjust_t.add_argument('-o',required=True,help='Directory for saving output')
    bulkDEGD_adjust_t.add_argument('-m',required=True,help='Group informations for analysis')
    bulkDEGD_adjust_t.add_argument('-s',required=True,help='Specify the species',choices=["pig"])
    bulkDEGD_adjust_t.add_argument('-p',help='Whether to draw plot or not',action='store_true')
    bulkDEGD_adjust_t.add_argument('-H',help='Data types for drawing heatmap plots',default='count')
    bulkDEGD_adjust_t.add_argument('-v',help='Data types for drawing volcano plots',default='qvalue')
    bulkDEGD_adjust_t.add_argument('-g',help='Interest gene for plot heatmap,split by ","',default='')
    bulkDEGD_adjust_t.add_argument('-t',help='Number of top differential expression genes for plot heatmap',default=30,type=int)
    bulkDEGD_adjust_t.add_argument('-fc',help='fc_threshold for vocano plot',default=1.5)

    # bulkRNA: 5. GeneOntologyAnalysis
    # 添加 bulkRNA GO 子命令
    bulkGO_parser = bulkRNA_subparsers.add_parser('enrichment', help='bulkRNA-seq QC')
    bulkGO_subparser = bulkGO_parser.add_subparsers(dest='GO_mode', required=True, help="Choose a mode to run the script")
    # enrichment_analysis_customer
    bulkGO_custom = bulkGO_subparser.add_parser('custom', help='Perform the enrichment analysis based on the customization gene sets')
    bulkGO_custom.add_argument('-i', required=True, help='path of a single column file containing genes of interests(SYMBOL)')
    bulkGO_custom.add_argument('-c', required=True, help='A customer file containing two columns, every row is consisted of two factors: the customer term name and gene of interests separated by tab')
    bulkGO_custom.add_argument('-o', required=True, help='Specify the absolute directory of output file, the ending does not require a /')
    # GO enrichment_analysis
    bulkGO = bulkGO_subparser.add_parser('GO', help='Perform the GO enrichment analysis based on the pubsihed GO/KEGG gene sets')
    bulkGO.add_argument('-i', required=True, help='path of a single column file containing genes of interests(SYMBOL)')
    bulkGO.add_argument('-s', required=True, help='Choose the functional annotation database of specific species: human, rat, mouse, chicken, or pig',
                        choices=["human", "rat","mouse","chicken","pig"])
    bulkGO.add_argument('-t', required=True, help='Define a threshold between 0~1, the GO terms will be merged into one once the gene sets matched to them exceed the threshold')
    bulkGO.add_argument('-o', required=True, help='Specify the absolute directory of output file, the ending does not require a /')
    # KEGG enrichment_analysis
    bulkKEGG = bulkGO_subparser.add_parser('KEGG', help='Perform the KEGG enrichment analysis based on the pubsihed GO/KEGG gene sets')
    bulkKEGG.add_argument('-i', required=True, help='path of a single column file containing genes of interests(SYMBOL)')
    bulkKEGG.add_argument('-s', required=True, help='Choose the functional annotation database of specific species: human, rat, mouse, chicken, or pig',
                        choices=["human", "rat","mouse","chicken","pig"])
    bulkKEGG.add_argument('-o', required=True, help='Specify the absolute directory of output file, the ending does not require a /')
    # plot_enrichment
    bulkGO_plot = bulkGO_subparser.add_parser('plot', help='Enrichment information of interests which can be')
    bulkGO_plot.add_argument('-i', required=True, help='Enrichment information of interests which can be')
    bulkGO_plot.add_argument('--color_high', required=True, help='Assign a color representing the maximum –lg(p.adjust) value')
    bulkGO_plot.add_argument('--color_low', required=True, help='Assign a color representing the minimum –lg(p.adjust) value')
    bulkGO_plot.add_argument('-o', required=True, help='Specify the absolute directory of output file, the ending does not require a /')

    # bulkRNA: 6. Similarity Analysis
    bulkSimilarity_parser = bulkRNA_subparsers.add_parser('similarity', help='bulkRNA-seq similarity analysis')
    # Required parameters
    bulkSimilarity_parser.add_argument('-b', required=True, help='Secify the base_input file (fpkm.csv), column name as sample, row name as gene')
    bulkSimilarity_parser.add_argument('-r', required=True, help='Specify the ref_input file (fpkm.csv), column name as sample, row name as gene')
    bulkSimilarity_parser.add_argument('-o', required=True, help='Specify the absolute directory of output file, the ending does not require a /')
    # Optional parameters
    bulkSimilarity_parser.add_argument('-x',default='NO',help='Specify the symbol representing repeated experiments of base data')
    bulkSimilarity_parser.add_argument('-y',default='NO',help='Specify the symbol representing repeated experiments of reference data')
    bulkSimilarity_parser.add_argument('-g',default='ALL',help='Specify the gene list (if not, fill in ALL)')

    # bulkRNA: 7. WGCNA
    bulkWGCNA_parser = bulkRNA_subparsers.add_parser('WGCNA', help='bulkRNA-seq WGCNA')
    bulkWGCNA_parser.add_argument('--count_file', required=True, help='The path to the file containing the gene expression count data.')
    bulkWGCNA_parser.add_argument('--pheno_file', required=True, help='A file path containing phenotypic information that describes the characteristics or experimental conditions of each sample')
    bulkWGCNA_parser.add_argument('--trait',required=True, help='Phenotypic information, calculating the correlation between gene modules and phenotypes')
    bulkWGCNA_parser.add_argument('-o',required=True,help='Output directory path')

    # bulkRNA: 8. AlternativeSplicingDetection
    bulkSplice_parser = bulkRNA_subparsers.add_parser('splice', help='bulkRNA-seq WGCNA')
    bulkSplice_subparser = bulkSplice_parser.add_subparsers(dest='Splice_mode', required=True, help="Choose a mode to run the script")
    # enrichment_analysis_customer.
    bulkSplice_enrich = bulkSplice_subparser.add_parser('detect', help='Detect the differentail splicing events by rMATs based on the bam file')
    bulkSplice_enrich.add_argument('--b1', required=True, help='Path to the b1 file containing BAM file paths separated by comma')
    bulkSplice_enrich.add_argument('--b2', required=True, help='Path to the b2 file containing BAM file paths separated by comma')
    bulkSplice_enrich.add_argument('--gtf', required=True, help='Path to the GTF file')
    bulkSplice_enrich.add_argument('-o', required=True, help='Output directory for the results')
    # draw_sashimiplot
    bulkSplice_draw = bulkSplice_subparser.add_parser('draw', help='Draw the sashimi plot to visualize the differential splicing events based on the bam files')
    bulkSplice_draw.add_argument('--b1', required=True, help='Path to the b1 file containing BAM file paths separated by comma')
    bulkSplice_draw.add_argument('--b2', required=True, help='Path to the b2 file containing BAM file paths separated by comma')
    bulkSplice_draw.add_argument('--gff', required=True, help='Path to the GFF file')
    bulkSplice_draw.add_argument('-o', required=True, help='Output directory for the results')
    bulkSplice_draw.add_argument('--interval', required=True, help='Interval of interests,chromosome:orientation:start:end. example: chr11:-:4752513:4854957')
    bulkSplice_draw.add_argument('--group', required=True, help='Path to the file containing group information of samples. example: WT: 1-3')

    # tools
    tools_parser = subparsers.add_parser('tools', help='some other useful tools')
    tools_subparser = tools_parser.add_subparsers(help='some other useful tools', title='tools', dest='tools_command')

    bw_parser = tools_subparser.add_parser('bw', help='Transfer bam file to bigwig file')
    bw_parser.add_argument('-i', required=True, help='Specify the absolute directory of input file (.bam), the ending does not require a /')
    bw_parser.add_argument('-o', required=True, help='Specify the absolute directory of output file, the ending does not require a /')
    bw_parser.add_argument('-p', help='Specify number of threads (default is 8)',default=8,type=int)

    te_parser = tools_subparser.add_parser('te', help='Quantify transposable element expression')
    te_parser.add_argument('-i', required=True, help='Specify the absolute directory of input file (.bam), the ending does not require a /')
    te_parser.add_argument('-o', required=True, help='Specify the absolute directory of output file, the ending does not require a /')
    te_parser.add_argument('-p', help='Specify number of threads (default is 8)',default=8,type=int)


    # scRNA
    # 创建 scRNA 子命令的子解析器
    scRNA_parser = subparsers.add_parser('scRNA', help='scRNA pipeline')
    scRNA_subparser = scRNA_parser.add_subparsers(help='scRNA pipeline', title='scRNA pipeline', dest='scRNA_command')

    # scRNA: 1. preprocess
    scRNApreprocess_parser = scRNA_subparser.add_parser('preprocess', help='scRNA preprocess command')
    scRNApreprocess_subparser = scRNApreprocess_parser.add_subparsers(dest='scRNA_preprocess_mode',required=True,help='Choose a mode')
    # cellranger
    scRNApreprocess_cellranger = scRNApreprocess_subparser.add_parser('cellranger',help='Perform cellranger preprocess')
    scRNApreprocess_cellranger.add_argument('-i', required=True, help='Input folder of fastq.gz files: X_1.fastq.gz,X_2.fastq.gz')
    scRNApreprocess_cellranger.add_argument('-s', required=True, help='Specify the species',choices=["pig",'test'])
    scRNApreprocess_cellranger.add_argument('-o', required=True, help='Output directory for the results')
    scRNApreprocess_cellranger.add_argument('-t',help='Specify number of threads (default is 8)',default=8,type=int)
    # DrSeq
    scRNApreprocess_drseq = scRNApreprocess_subparser.add_parser('drseq',help='Perform DrSeq preprocess')
    scRNApreprocess_drseq.add_argument('-i', required=True, help='Input folder of fastq files: X_1.fastq,X_2.fastq')
    scRNApreprocess_drseq.add_argument('-s', required=True, help='Specify the species',choices=["pig"])
    scRNApreprocess_drseq.add_argument('-o', required=True, help='Output directory for the results')
    scRNApreprocess_drseq.add_argument('-t',help='Specify number of threads (default is 8)',default=8,type=int)
    # smartseq
    scRNApreprocess_smartseq = scRNApreprocess_subparser.add_parser('smartseq',help='Perform smartseq preprocess')
    scRNApreprocess_smartseq.add_argument('-i', required=True, help='Specify the absolute directory of input file (fastq.gz), the ending does not require a /')
    scRNApreprocess_smartseq.add_argument('-o', required=True, help='Specify the absolute directory of output file, the ending does not require a /')
    scRNApreprocess_smartseq.add_argument('-s', required=True, help='Specify the species',
                                      choices=["mm9", "mm10","pig"])
    scRNApreprocess_smartseq.add_argument('-bw', action='store_true', help='If specified, perform a bam -> bigwig file transfer')
    scRNApreprocess_smartseq.add_argument('-te', action='store_true', help='If specified, perform a TE expression analysis')
    scRNApreprocess_smartseq.add_argument('-t', help='Specify number of threads (default is 8)',default=8,type=int)
    scRNApreprocess_smartseq.add_argument('-l', help='Single sequencing, 2:Pair sequencing (default is 2)',default=2)
    # singleron
    scRNApreprocess_singleron = scRNApreprocess_subparser.add_parser('singleron',help='Perform singleron preprocess')
    scRNApreprocess_singleron.add_argument('-i', required=True, help='Specify the input fastq file directory, fastqs file format: X_1.fastq.gz, X_2.fastq.gz')
    scRNApreprocess_singleron.add_argument('-o', required=True, help='Specify the absolute directory of output file, the ending does not require a /')
    scRNApreprocess_singleron.add_argument('-s', required=True, help='Specify the species',
                                      choices=["galGal"])
    scRNApreprocess_singleron.add_argument('-t', help='Specify number of threads (default is 8)',default=8,type=int)

    # scRNA: 3. clustering
    scRNAclustering_parser = scRNA_subparser.add_parser('cluster', help='scRNA clustering command')
    scRNAclustering_parser.add_argument('-i',required=True, help='Input files directory')
    scRNAclustering_parser.add_argument('-o',required=True, help='Ouput directory')
    scRNAclustering_parser.add_argument('-n',required=True, help='Output prefix name')

    scRNAclustering_parser.add_argument('--min_cells', type=int, default=3, help="clustering minimum cells, default: 3")
    scRNAclustering_parser.add_argument('--min_features', type=int, default=200, help="clustering maximum cells, default: 200")
    scRNAclustering_parser.add_argument('--nFeature_RNA_min', type=int, default=200, help="clustering feature minimum, default: 200")
    scRNAclustering_parser.add_argument('--nFeature_RNA_max', type=int, default=2500, help="clustering feature maximum, default: 2500")
    scRNAclustering_parser.add_argument('--mt_percent', type=float, default=5, help="clustering mitochondria percent threshold, default: 5")
    scRNAclustering_parser.add_argument('--variable_features', type=int, default=2000, help="clustering variable features, default: 2000")
    scRNAclustering_parser.add_argument('--dim_prefer', type=int, default=30, help="clustering dimesions, default: 30")
    scRNAclustering_parser.add_argument('--res', type=float, default=1, help="clustering resolution, default: 1")

    # scRNA: 7. trajectory
    scRNAtrajectory_parser = scRNA_subparser.add_parser('trajectory', help='scRNA trajectory command')
    scRNAtrajecory_subparser = scRNAtrajectory_parser.add_subparsers(dest='scRNA_trajectory_mode',required=True,help='Choose a mode')
    # monocle3
    scRNAtrajecory_monocle = scRNAtrajecory_subparser.add_parser('monocle',help='Perform monocle trajectory analysis')
    scRNAtrajecory_monocle.add_argument('-i', required=True, help='Input rds file, clustered and annotated seurat object')
    scRNAtrajecory_monocle.add_argument('-m', required=True, help='Meta data for trajectory analysis, e.g., celltype')
    scRNAtrajecory_monocle.add_argument('-o', required=True, help='Output directory for the results')
    scRNAtrajecory_monocle.add_argument('-s', required=True, help='Start point for trajectory analysis, e.g., NK')
    # DiffusionMap
    scRNAtrajecory_diffusionmap = scRNAtrajecory_subparser.add_parser('diffusion',help='Perform diffusion map trajectory analysis')
    scRNAtrajecory_diffusionmap.add_argument('-i', required=True, help='Input rds file, clustered and annotated seurat object')
    scRNAtrajecory_diffusionmap.add_argument('-o', required=True, help='Output directory for the results')
    scRNAtrajecory_diffusionmap.add_argument('-m', required=True, help='Meta data for trajectory analysis, e.g., celltype')
    scRNAtrajecory_diffusionmap.add_argument('-s', required=True, help='Start point for trajectory analysis, e.g., NK')
    # velocity
    scRNAtrajecory_velocity = scRNAtrajecory_subparser.add_parser('velocity',help='Perform velocity trajectory analysis')
    scRNAtrajecory_velocity.add_argument('--cellranger', required=True, help='Cellranger output directory')
    scRNAtrajecory_velocity.add_argument('--gtf', required=True, help='GTF file for velocity analysis')
    scRNAtrajecory_velocity.add_argument('--rds', required=True, help='rds file, clustered and annotated seurat object for velocity analysis')
    scRNAtrajecory_velocity.add_argument('--genes', type=str,default='',help='Interested genes, separate by ","')
    scRNAtrajecory_velocity.add_argument('-o', required=True, help='Output directory for the results')
    # CStreet
    scRNAtrajecory_cstreet = scRNAtrajecory_subparser.add_parser('cstreet',help='Perform CStreet trajectory analysis')
    scRNAtrajecory_cstreet.add_argument("-i", type=str, help="Input expression matrixes, separate by ',' ")
    scRNAtrajecory_cstreet.add_argument("-s", type=str, help="Input cell state file, separate by ',' ")
    scRNAtrajecory_cstreet.add_argument("-n", type=str,default='',help="Output directory name")
    scRNAtrajecory_cstreet.add_argument("-o", type=str, help="Output directory")

    # 解析命令行参数
    args = parser.parse_args()

    # 判断并执行 bulkRNA 相关子命令
    if args.omics == 'bulkRNA':
        if args.bulkRNA_command == 'preprocess':
            bulkRNA_preprocess(args)
        elif args.bulkRNA_command == 'QC':
            bulkRNA_QC(args)
        elif args.bulkRNA_command == 'DEG':
            if args.DEG_mode is None:
                print("\nError: Must specify a mode ('deseq2' , 'gfold' or 'adjusted_t').")
            bulkRNA_DEG(args)
        elif args.bulkRNA_command == 'enrichment':
            bulkRNA_GO(args)
        elif args.bulkRNA_command == 'similarity':
            bulkRNA_Similarity(args)
        elif args.bulkRNA_command == 'WGCNA':
            bulkRNA_WGCNA(args)
        elif args.bulkRNA_command == 'splice':
            bulkRNA_AlternativeSpliceDetection(args)
            
    if args.omics == 'tools':
        if args.tools_command == 'bw':
            bw(args)
        elif args.tools_command == 'te':
            te(args)

    if args.omics == 'scRNA':
        if args.scRNA_command == 'preprocess':
            scRNA_preprocess(args)
        elif args.scRNA_command == 'cluster':
            scRNA_cluster(args)
        elif args.scRNA_command == 'trajectory':
            scRNA_trajectory(args)

if __name__ == "__main__":
    main()



