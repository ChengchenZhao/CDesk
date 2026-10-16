"""
ATAC-seq params
"""

def register_subparser(subparsers, config):
    ATAC_parser = subparsers.add_parser('ATAC', help='ATAC pipeline')
    ATAC_subparsers = ATAC_parser.add_subparsers(
        help='ATAC pipeline', title='ATAC pipeline', dest='ATAC_command'
    )

    # ATAC: 1. preprocess
    ATACpreprocess_parser = ATAC_subparsers.add_parser('preprocess', help='ATAC preprocess command')
    ATACpreprocess_parser.add_argument('-i','--input', required=True, help='Input fq information file')
    ATACpreprocess_parser.add_argument('-o','--output', required=True, help='Output directory')
    ATACpreprocess_parser.add_argument('-s','--species', required=True, help='Specify species', choices=list(config['data'].keys()))
    ATACpreprocess_parser.add_argument('-t','--thread', type=int, default=8, help='Number of threads, default: 8')

    # ATAC: 2. cor
    ATAC_QC_parser = ATAC_subparsers.add_parser('correlation', help='ATAC sample correlation')
    ATAC_QC_parser.add_argument("-i","--input", required=True, type=str, help="The input ATAC sample information file")
    ATAC_QC_parser.add_argument("-s","--species", required=True, type=str, help="Species", choices=list(config['data'].keys()))
    ATAC_QC_parser.add_argument("--bin", type=int, default=100000, help="Bin size for genomic binning, default: 100000")
    ATAC_QC_parser.add_argument("--step", type=int, default=1, help="Peak file partition step, default: 1")
    ATAC_QC_parser.add_argument("-o","--output", required=True, type=str, help="Output directory")
    ATAC_QC_parser.add_argument("-t","--thread", type=int, default=30, help="Threads to run, default: 30")
    ATAC_QC_parser.add_argument('--width', default=12, type=int, help='The plot width, default: 12')
    ATAC_QC_parser.add_argument('--height', default=8, type=int, help='If plot height, default:8')

    # ATAC: 3. Dynamic
    ATAC_dynamic_parser = ATAC_subparsers.add_parser('dynamic', help='ATAC dynamic pattern analysis')
    ATAC_dynamic_parser.add_argument("-i","--input", required=True, type=str, help="The input sample information file")
    ATAC_dynamic_parser.add_argument("--remove_cluster", default='', type=str, help="The cluster to remove, default: ''")
    ATAC_dynamic_parser.add_argument("-o","--output", required=True, type=str, help="The output directory")
    ATAC_dynamic_parser.add_argument("--important_index", default=0, type=int, help="Specify the group to align, 1: the latter one in meta file, 0(default): the former one in meta file")
    ATAC_dynamic_parser.add_argument("-s","--species", required=True, type=str, choices=list(config['data'].keys()), help="Species")

    # ATAC: 4. motif
    ATAC_motif_parser = ATAC_subparsers.add_parser('motif', help='ATAC motif analysis')
    ATAC_motif_subparser = ATAC_motif_parser.add_subparsers(dest='motif_mode', required=True, help="Choose to analyze or plot")
    # homer find
    ATAC_Find_homer_subparser = ATAC_motif_subparser.add_parser('homer', help='Perform homer motif finding analysis')
    ATAC_Find_homer_subparser.add_argument("-i","--input", required=True, type=str, help="The input gene/peak information file")
    ATAC_Find_homer_subparser.add_argument("-o","--output", required=True, type=str, help="Output directory")
    ATAC_Find_homer_subparser.add_argument("-s","--species", required=True, help='Available gene in HOMER or provide the path to genome FASTA files for peak mode', type=str)
    ATAC_Find_homer_subparser.add_argument("--mode", required=True, type=str, help="Analysis mode: peak/gene", choices=['gene','peak'])
    ATAC_Find_homer_subparser.add_argument("-t","--thread", default=20, type=int, help="Number of threads, default: 20")
    ATAC_Find_homer_subparser.add_argument("--start", default=-300, type=int, help="offset from TSS, default=-300")
    ATAC_Find_homer_subparser.add_argument("--end", default=50, type=int, help="offset from TSS, default=50")
    # meme find
    ATAC_Find_meme_subparser = ATAC_motif_subparser.add_parser('meme', help='Perform meme motif finding analysis')
    ATAC_Find_meme_subparser.add_argument("-i","--input", required=True, help="The fasta input information file")
    ATAC_Find_meme_subparser.add_argument("-o","--output", required=True, type=str, help="Output directory")
    ATAC_Find_meme_subparser.add_argument("-t","--thread", default=10, type=int, help="Number of threads, default: 10")
    ATAC_Find_meme_subparser.add_argument("--nmotifs", default=3, type=int, help="Number of motifs to extract, default: 3")
    ATAC_Find_meme_subparser.add_argument("-db","--database", required=True, type=str, help="Database file")
    # homer heatmap
    ATAC_homer_heatmap_subparser = ATAC_motif_subparser.add_parser('homer_heatmap', help='Plot the homer finding result heatmap')
    ATAC_homer_heatmap_subparser.add_argument("-i","--input", required=True, type=str, help="The homer results information csv file")
    ATAC_homer_heatmap_subparser.add_argument("-o","--output", required=True, type=str, help="The output directory")
    ATAC_homer_heatmap_subparser.add_argument("--motif", default='no', type=str, help="The specified motifs txt file")
    ATAC_homer_heatmap_subparser.add_argument("--width", default=10, type=float, help="The plot width, default: 10")
    ATAC_homer_heatmap_subparser.add_argument("--height", default=8, type=float, help="The plot height, default: 8")
    # enrich heatmap
    ATAC_enrich_heatmap_subparser = ATAC_motif_subparser.add_parser('enrich_heatmap', help='Perform the motif enriched heatmap')
    ATAC_enrich_heatmap_subparser.add_argument("--motif_bed", required=True, type=str, help="The motif bed file information csv file")
    ATAC_enrich_heatmap_subparser.add_argument("--bed", required=True, type=str, help="The input bed file information csv file")
    ATAC_enrich_heatmap_subparser.add_argument("-o","--output", required=True, type=str, help="The output directory")
    ATAC_enrich_heatmap_subparser.add_argument("-t","--thread", default=12, type=int, help="The number of threads")
    ATAC_enrich_heatmap_subparser.add_argument("-s","--species", required=True, type=str, help="Species")

    # ATAC: 5. Accessbility
    ATAC_accessbility_parser = ATAC_subparsers.add_parser('accessbility', help='ATAC accessbility analysis')
    ATAC_accessbility_parser.add_argument("--bw", required=True, help="The input bw txt file")
    ATAC_accessbility_parser.add_argument("-o","--output", required=True, help="Output directory")
    ATAC_accessbility_parser.add_argument("--mode", required=True, choices=["region", "gene"], help="Analysis mode: region/gene")
    ATAC_accessbility_parser.add_argument("-i","--input", required=True, help="Input gene region or gene list txt file")
    ATAC_accessbility_parser.add_argument("--species", help="Species", choices=list(config['data'].keys()))
    ATAC_accessbility_parser.add_argument("--method", choices=["sample", "region_list"], default="sample", help="Boxplot plotting method：sample/region list")
    ATAC_accessbility_parser.add_argument("--center", required=True, help="Average profile alignment：region mode:peak/position，gene mode:TSS/TTS/whole")
    ATAC_accessbility_parser.add_argument("--promoter_size", type=int, default=5000, help="Gene promoter region size, default: TSS upstream and downstream 5kb")
    ATAC_accessbility_parser.add_argument("--region_list", default='', help="Region list file for box plot")

    # ATAC: 6. IGV snapshot
    ATAC_igv_parser = ATAC_subparsers.add_parser('igv', help='ATAC IGV snapshot')
    ATAC_igv_parser.add_argument("-i","--input", required=True, type=str, help="The input list file")
    ATAC_igv_parser.add_argument("-o","--output", required=True, type=str, help="Output directory")
    ATAC_igv_parser.add_argument("--genome", required=True, type=str, help="Select the genome id")
    ATAC_igv_parser.add_argument("--display_mode", default='collapse', help="Display mode", choices=["collapse","expand","squish"])
    ATAC_igv_parser.add_argument("--region", required=True, help="Specify the regions, separate by ,")
    ATAC_igv_parser.add_argument("--type", default='png', type=str, choices=['svg','png'], help="svg or png output")

    # ATAC: 7. Similarity
    ATAC_correlation_parser = ATAC_subparsers.add_parser('similarity', help='ATAC sample similarity analysis')
    ATAC_correlation_parser.add_argument("--bw1", required=True, help='Path to the first bw file')
    ATAC_correlation_parser.add_argument("--bw2", required=True, help='Path to the second bw file')
    ATAC_correlation_parser.add_argument("-o","--output", required=True, help='Output directory')
    ATAC_correlation_parser.add_argument("-t","--thread", type=int, default=30, help='Number of threads')
    ATAC_correlation_parser.add_argument("--region", help='Optional BED file for specific regions')

    return ATAC_parser
