"""
bulkRNA functions
"""

def register_subparser(subparsers, config):
    bulkRNA_parser = subparsers.add_parser('bulkRNA', help='bulkRNA-seq pipeline')
    bulkRNA_subparsers = bulkRNA_parser.add_subparsers(
        help='bulkRNA-seq pipeline', title='bulkRNA-seq pipeline', dest='bulkRNA_command'
    )

    # bulkRNA: 1. preprocess
    bulkPreprocessing_parser = bulkRNA_subparsers.add_parser('preprocess', help='bulkRNA-seq preprocessing')
    bulkPreprocessing_parser.add_argument('--input','-i', required=True, help='The input sequencing data information file')
    bulkPreprocessing_parser.add_argument('--output','-o', required=True, help='Specify the output directory')
    bulkPreprocessing_parser.add_argument('--species','-s', required=True, help='Specify the species',
                                          choices=list(config['data'].keys()))
    bulkPreprocessing_parser.add_argument('-bw', action='store_true', help='If specified, perform a bam -> bigwig file transfer (optional)')
    bulkPreprocessing_parser.add_argument('-te', action='store_true', help='If specified, perform a TE expression analysis (optional)')
    bulkPreprocessing_parser.add_argument('--thread','-t', help='Specify number of threads (default is 8)', default=8, type=int)

    # bulkRNA: 2. cor
    bulkQC_parser = bulkRNA_subparsers.add_parser('correlation', help='bulkRNA-seq sample correlation')
    bulkQC_parser.add_argument('--input','-i', required=True, help='Specify the input gene expression file (.csv), column name as sample, row name as gene')
    bulkQC_parser.add_argument('--output','-o', required=True, help='Specify the output directory')
    bulkQC_parser.add_argument('--group', required=True, help='The grouping file(csv,txt,tsv)')
    bulkQC_parser.add_argument('--batch', default='no', help='Specify the batch effect removing method, default: no',
                               choices=['no','removeBatchEffect','combat'])
    bulkQC_parser.add_argument('--width', default=8, type=int, help='The plotting width, default: 8')
    bulkQC_parser.add_argument('--height', default=6, type=int, help='The plotting height, default: 6')

    # bulkRNA: 3. DEG
    bulkDEGD_parser = bulkRNA_subparsers.add_parser('DEG', help='bulkRNA-seq differential gene expression analysis')
    bulkDEGD_subparser = bulkDEGD_parser.add_subparsers(dest='DEG_mode', required=True, help="Choose a mode to run the script")
    # DEseq2
    bulkDEGD_deseq2 = bulkDEGD_subparser.add_parser('deseq2', help='Perform DEseq2 differential gene expression analysis')
    bulkDEGD_deseq2.add_argument('-i','--input', required=True, help='Input count matrix for DESeq2')
    bulkDEGD_deseq2.add_argument('-o','--output', required=True, help='Output directory')
    bulkDEGD_deseq2.add_argument('--group', required=True, help='Grouping file')
    bulkDEGD_deseq2.add_argument('-p', help='Whether to draw plot or not', action='store_true')
    bulkDEGD_deseq2.add_argument('--gene', help='Interested genes file', default='NO')
    bulkDEGD_deseq2.add_argument('--top_num', help='Number of top differential expression genes for heatmap if no gene file provided, default: 500', default=500, type=int)
    bulkDEGD_deseq2.add_argument('--fc', help='fold change threshold of DEGs, default: 1', type=float, default=1)
    bulkDEGD_deseq2.add_argument('--pval', help='p_adjusted threshold of DEGs, default: 0.05', type=float, default=0.05)
    bulkDEGD_deseq2.add_argument('--width', help='Plot width, default: 5', type=float, default=5)
    bulkDEGD_deseq2.add_argument('--height', help='Plot height, default: 5', type=float, default=5)
    # adjusted_t
    bulkDEGD_adjust_t = bulkDEGD_subparser.add_parser('adjusted_t', help='Perform adjusted_t differential gene expression analysis')
    bulkDEGD_adjust_t.add_argument('-i','--input', required=True, help='Input matrix for adjusted_t')
    bulkDEGD_adjust_t.add_argument('-o','--output', required=True, help='Output directory')
    bulkDEGD_adjust_t.add_argument('--group', required=True, help='Grouping file')
    bulkDEGD_adjust_t.add_argument('-p', help='Whether to draw plot or not', action='store_true')
    bulkDEGD_adjust_t.add_argument('--gene', help='Interested genes file', default='NO')
    bulkDEGD_adjust_t.add_argument('--top_num', help='Number of top differential expression genes for heatmap if no gene file provided, default: 500', default=500, type=int)
    bulkDEGD_adjust_t.add_argument('--fc', help='fold change threshold of DEGs, defaullt: 1', type=float, default=1)
    bulkDEGD_adjust_t.add_argument('--pval', help='p_adjusted threshold of DEGs, default: 0.05', type=float, default=0.05)
    bulkDEGD_adjust_t.add_argument('--width', help='Plot width, default: 5', type=float, default=5)
    bulkDEGD_adjust_t.add_argument('--height', help='Plot height, default: 5', type=float, default=5)
    # gfold
    bulkDEGD_gfold = bulkDEGD_subparser.add_parser('gfold', help='Perform adjusted_t differential gene expression analysis')
    bulkDEGD_gfold.add_argument('-i','--input', required=True, help='Input bam grouping file')
    bulkDEGD_gfold.add_argument('-o','--output', required=True, help='Output directory')
    bulkDEGD_gfold.add_argument('-s','--species', required=True, help='Specify the species', choices=list(config['data'].keys()))
    bulkDEGD_gfold.add_argument('-p', help='Whether to draw plot or not', action='store_true')
    bulkDEGD_gfold.add_argument('-t','--thread', help='Number of threads, default: 20', type=int, default=20)
    bulkDEGD_gfold.add_argument('--gene', help='Interested gene file', default='NO')
    bulkDEGD_gfold.add_argument('--top_num', help='Number of top differential expression genes for heatmap if no gene file provided, default: 500', default=500, type=int)
    bulkDEGD_gfold.add_argument('--gfold', help='GFold threshold of DEGs, default: 1', type=float, default=1)
    bulkDEGD_gfold.add_argument('--width', help='Plot width, default: 5', type=float, default=5)
    bulkDEGD_gfold.add_argument('--height', help='Plot height, default: 5', type=float, default=5)

    # bulkRNA: 4. enrich
    bulkEnrich_parser = bulkRNA_subparsers.add_parser('enrich', help='bulkRNA-seq gene functional enrichment analysis')
    bulkEnrich_subparser = bulkEnrich_parser.add_subparsers(dest='enrich_mode', required=True, help="Choose a mode to run the script")
    # enrich analyze
    bulk_enrich = bulkEnrich_subparser.add_parser('analyze', help='Perform the functional enrichment analysis')
    bulk_enrich.add_argument('-i', '--input', required=True, help='The single column file of interested genes (SYMBOL)')
    bulk_enrich.add_argument('--type', help='Specify the analysis type, single/multiple samples', required=True, choices=["single","multi"])
    bulk_enrich.add_argument('-s','--species', help='Species', default='no', choices=["human", "rat","mouse","chicken","pig","tair","no"])
    bulk_enrich.add_argument('-o','--output', required=True, help='Output directory')
    bulk_enrich.add_argument('--custom', default='no', help='The reference customer file containing two columns for custom analysis (every row is consisted of two factors: the customer term name and gene of interests separated by tab)')
    bulk_enrich.add_argument('--width', default='10', type=int, help='Plot width, default: 10')
    bulk_enrich.add_argument('--height', default='6', type=int, help='Plot height, default: 6')
    # enrich plot
    bulk_plot = bulkEnrich_subparser.add_parser('plot', help='Plot the functional enrichment result')
    bulk_plot.add_argument('-i','--input', required=True, help='The enrichment result file')
    bulk_plot.add_argument('-o','--output', required=True, help='Output directory')
    bulk_plot.add_argument('--type', help='Specify the analysis type, single/multiple samples', required=True, choices=["single","multi"])
    bulk_plot.add_argument('--width', default='10', type=int, help='Plot width, default: 10')
    bulk_plot.add_argument('--height', default='6', type=int, help='Plot height, default: 6')
    # GSEA
    bulkGSEA = bulkEnrich_subparser.add_parser('GSEA', help='bulkRNA-seq GSEA analysis')
    bulkGSEA.add_argument("-i","--input", required=True, type=str, help="The ranked gene list file (e.g., DEG result file)")
    bulkGSEA.add_argument("-o","--output", required=True, type=str, help="Output directory")
    bulkGSEA.add_argument("-s","--species", required=True, type=str, help="Specify the species", choices=['human','mouse','pig','chicken','rat','tair'])
    bulkGSEA.add_argument("--path", type=str, default='NO', help="Specify the paths to plot")
    bulkGSEA.add_argument("--cols", default='gene_name,log2FoldChange', type=str, help="Columns used, default:gene_name,log2FoldChange")
    bulkGSEA.add_argument('--width', default=7, type=float, help='Plot width, default: 7')
    bulkGSEA.add_argument('--height', default=7, type=float, help='Plot height, default: 7')

    # bulkRNA: 5. Similarity Analysis
    bulkSimilarity_parser = bulkRNA_subparsers.add_parser('similarity', help='bulkRNA-seq similarity analysis')
    bulkSimilarity_parser.add_argument('-i','--input', required=True, help='The gene expression list file')
    bulkSimilarity_parser.add_argument('-o','--output', required=True, help='Output directory')
    bulkSimilarity_parser.add_argument('--group', required=True, help='The grouping file')
    bulkSimilarity_parser.add_argument('--gene', default='All', help='Specify the gene list, default: ALL')
    bulkSimilarity_parser.add_argument('--batch', default='no', help='Specify the batch effect removing method', choices=['no','removeBatchEffect','combat'])
    bulkSimilarity_parser.add_argument('--width', default=10, type=int, help='The plotting width, default: 10')
    bulkSimilarity_parser.add_argument('--height', default=8, type=int, help='If plotting height, default:8')

    # bulkRNA: 6. Clustering
    bulkClustering_parser = bulkRNA_subparsers.add_parser('cluster', help='bulkRNA-seq gene clustering analysis')
    bulkClustering_subparser = bulkClustering_parser.add_subparsers(dest='cluster_mode', required=True, help="Choose to use Kmean or WGCNA")
    # WGCNA
    bulkWGCNA_parser = bulkClustering_subparser.add_parser('WGCNA', help='bulkRNA-seq WGCNA analysis')
    bulkWGCNA_parser.add_argument('-i','--input', required=True, help='The input gene expression matrix file.')
    bulkWGCNA_parser.add_argument('--pheno', required=True, help='The sample trait information file.')
    bulkWGCNA_parser.add_argument('--trait', required=True, help='The phenotypes specified to calculate correlations with gene modules.')
    bulkWGCNA_parser.add_argument('-o','--output', required=True, help='Output directory')
    # Kmean
    bulkKmean_parser = bulkClustering_subparser.add_parser('kmean', help='bulkRNA-seq Knmean gene clustering')
    bulkKmean_parser.add_argument('-i','--input', required=True, help='The input gene expression matrix file.')
    bulkKmean_parser.add_argument('--group', required=True, help='The grouping csv file')
    bulkKmean_parser.add_argument('--cluster', type=int, default=6, help='Number of clusters for K-means clustering, default: 6')
    bulkKmean_parser.add_argument('-o','--output', required=True, help='Output directory')
    bulkKmean_parser.add_argument('--method', default='correlation', help='Distance measure to be used for clustering, default: correlation',
                                  choices=["euclidean", "maximum", "manhattan", "canberra", "binary", "pearson", "abspearson", "abscorrelation", "correlation","spearman", "kendall"])
    bulkKmean_parser.add_argument('--height', default=10, help='Output image height, default: 10', type=float)
    bulkKmean_parser.add_argument('--width', default=8, help='Output image width, defaullt: 8', type=float)
    bulkKmean_parser.add_argument('--gene', default='no', help='Spcefify the gene list txt file (optional)')

    # bulkRNA: 7. AlternativeSplicingDetection
    bulkSplice_parser = bulkRNA_subparsers.add_parser('splice', help='bulkRNA-seq alternative splicing detection')
    bulkSplice_subparser = bulkSplice_parser.add_subparsers(dest='Splice_mode', required=True, help="Choose a mode to run the script")
    # Splice detect
    bulkSplice_detect = bulkSplice_subparser.add_parser('detect', help='Detect the differentail splicing events by rMATs based on the bam file')
    bulkSplice_detect.add_argument('--b1', type=str, required=True, help="The first BAM files txt")
    bulkSplice_detect.add_argument('-s','--species', type=str, required=True, choices=list(config['data'].keys()), help="Specify the species")
    bulkSplice_detect.add_argument('-o','--output', type=str, required=True, help="The output directory")
    bulkSplice_detect.add_argument('--b2', type=str, default='No', help="The second BAM files txt")
    bulkSplice_detect.add_argument('--type', choices=['paired','single'], default='paired', help="Type of read: paired/single, default:paired")
    bulkSplice_detect.add_argument('-t','--thread', type=int, default=10, help="Number of threads, default: 10")
    bulkSplice_detect.add_argument('--length', type=int, default=150, help="The length of each read, default: 150")
    bulkSplice_detect.add_argument('--variable_read_length', action='store_true', help="Allow reads with lengths that differ from read Length")
    bulkSplice_detect.add_argument('--allow_clipping', action='store_true', help="Allow alignments with soft or hard clipping to be used")
    # Splice draw
    bulkSplice_draw = bulkSplice_subparser.add_parser('draw', help='Draw the sashimi plot to visualize the differential splicing events based on the bam files')
    bulkSplice_draw.add_argument('--b1', type=str, required=True, help="The first BAM files txt")
    bulkSplice_draw.add_argument('-s','--species', type=str, required=True, help="Specify the species")
    bulkSplice_draw.add_argument('-o','--output', type=str, required=True, help="The output directory")
    bulkSplice_draw.add_argument('--region', type=str, required=True, help="The genome region coordinates: format:{chromosome}:{strand}:{start}:{end}")
    bulkSplice_draw.add_argument('--l1', type=str, help="The first sample name", default='sample1')
    bulkSplice_draw.add_argument('--l2', type=str, help="The second sample name", default='sample2')
    bulkSplice_draw.add_argument('--group', type=str, default='No', help="The path to a .gf file which groups the replicates")
    bulkSplice_draw.add_argument('--b2', type=str, default='No', help="The second BAM files txt")
    bulkSplice_draw.add_argument('--exon_s', type=int, default=1, help="How much to scale down exon, default: 1")
    bulkSplice_draw.add_argument('--intron_s', type=int, default=1, help="How much to scale down introns, default: 1")

    return bulkRNA_parser
