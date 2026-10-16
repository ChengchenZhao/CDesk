"""
ChIPseqCUTTag params
"""

def register_subparser(subparsers, config):
    ChIPseq_CUTTag_parser = subparsers.add_parser('ChIPseqCUTTag', help='ChIPseq&CUTTag pipeline')
    ChIPseq_CUTTag_subparsers = ChIPseq_CUTTag_parser.add_subparsers(
        help='ChIPseq&CUTTag pipeline', title='ChIPseq&CUTTag pipeline', dest='ChIPseq_CUTTag_command'
    )

    # ChIPseq&CUTTag: 1. Preprocess
    ChIPseq_CUTTag_preprocess_parser = ChIPseq_CUTTag_subparsers.add_parser('preprocess', help='ChIPseq&CUTTag preprocess command')
    ChIPseq_CUTTag_preprocess_parser.add_argument('-i','--input', required=True, help='Input fq information file')
    ChIPseq_CUTTag_preprocess_parser.add_argument('-o','--output', required=True, help='Output directory')
    ChIPseq_CUTTag_preprocess_parser.add_argument('-s','--species', required=True, help='Specify species', choices=list(config['data'].keys()))
    ChIPseq_CUTTag_preprocess_parser.add_argument('-t','--thread', type=int, default=8, help='Number of threads, default: 8')

    # ChIPseq_CUTTag: 2. QC (Correlation)
    ChIPseq_CUTTag_QC_parser = ChIPseq_CUTTag_subparsers.add_parser('correlation', help='ChIPseq_CUTTag sample correlation')
    ChIPseq_CUTTag_QC_parser.add_argument("-i","--input", required=True, type=str, help="The input ChIPseq_CUTTag sample information file")
    ChIPseq_CUTTag_QC_parser.add_argument("-s","--species", required=True, type=str, help="Species", choices=list(config['data'].keys()))
    ChIPseq_CUTTag_QC_parser.add_argument("--bin", type=int, default=100000, help="Bin size for genomic binning, default: 100000")
    ChIPseq_CUTTag_QC_parser.add_argument("--step", type=int, default=1, help="Peak file partition step, default: 1")
    ChIPseq_CUTTag_QC_parser.add_argument("-o","--output", required=True, type=str, help="Output directory")
    ChIPseq_CUTTag_QC_parser.add_argument("-t","--thread", type=int, default=30, help="Threads to run, default: 30")
    ChIPseq_CUTTag_QC_parser.add_argument('--width', default=12, type=int, help='The plot width, default: 12')
    ChIPseq_CUTTag_QC_parser.add_argument('--height', default=8, type=int, help='If plot height, default:8')

    # ChIPseq&CUTTag: 3. Pattern
    ChIPseq_CUTTag_pattern_parser = ChIPseq_CUTTag_subparsers.add_parser('pattern', help='ChIPseq&CUTTag pattern analysis')
    ChIPseq_CUTTag_pattern_parser.add_argument("-bw", required=True, type=str, help="The input bw files list file")
    ChIPseq_CUTTag_pattern_parser.add_argument("-bed", required=True, type=str, help="The bed genomic regions to be analyzed list file")
    ChIPseq_CUTTag_pattern_parser.add_argument("-o","--output", required=True, type=str, help="Output directory")
    ChIPseq_CUTTag_pattern_parser.add_argument("--mode", default='reference-point', help="Analysis mode, default: reference-point", choices=["reference-point","scale-regions"])
    ChIPseq_CUTTag_pattern_parser.add_argument("--height", help="Heatmap height, default: 10", type=int, default=10)
    ChIPseq_CUTTag_pattern_parser.add_argument("-a","--upstream", help="Length to extend upstream of the genomic region start site (bp), default: 2000", type=int, default=2000)
    ChIPseq_CUTTag_pattern_parser.add_argument("-b","--downstream", help="Length to extend downstream of the genomic region end site (bp), default: 2000", type=int, default=2000)
    ChIPseq_CUTTag_pattern_parser.add_argument("-l","--length", help="Scaled genomic region length (bp), default: 3000", type=int, default=3000)
    ChIPseq_CUTTag_pattern_parser.add_argument("-t","--thread", help="Number of processors, default: 10", type=int, default=10)
    ChIPseq_CUTTag_pattern_parser.add_argument("--region", help="The reference point for reference-point, default: center", default='center', choices=['TSS','TES','center'])

    # ChIPseq_CUTTag: 4. motif
    ChIPseq_CUTTag_motif_parser = ChIPseq_CUTTag_subparsers.add_parser('motif', help='ChIPseq_CUTTag motif analysis')
    ChIPseq_CUTTag_motif_subparser = ChIPseq_CUTTag_motif_parser.add_subparsers(dest='motif_mode', required=True, help="Choose to analyze or plot")
    # homer find
    ChIPseq_CUTTag_Find_homer_subparser = ChIPseq_CUTTag_motif_subparser.add_parser('homer', help='Perform homer motif finding analysis')
    ChIPseq_CUTTag_Find_homer_subparser.add_argument("-i","--input", required=True, type=str, help="The input gene/peak information file")
    ChIPseq_CUTTag_Find_homer_subparser.add_argument("-o","--output", required=True, type=str, help="Output directory")
    ChIPseq_CUTTag_Find_homer_subparser.add_argument("-s","--species", required=True, help='Available gene in HOMER or provide the path to genome FASTA files for peak mode', type=str)
    ChIPseq_CUTTag_Find_homer_subparser.add_argument("--mode", required=True, type=str, help="Analysis mode: peak/gene", choices=['gene','peak'])
    ChIPseq_CUTTag_Find_homer_subparser.add_argument("-t","--thread", default=20, type=int, help="Number of threads, default: 20")
    ChIPseq_CUTTag_Find_homer_subparser.add_argument("--start", default=-300, type=int, help="offset from TSS, default=-300")
    ChIPseq_CUTTag_Find_homer_subparser.add_argument("--end", default=50, type=int, help="offset from TSS, default=50")
    # meme find
    ChIPseq_CUTTag_Find_meme_subparser = ChIPseq_CUTTag_motif_subparser.add_parser('meme', help='Perform meme motif finding analysis')
    ChIPseq_CUTTag_Find_meme_subparser.add_argument("-i","--input", required=True, help="The fasta input information file")
    ChIPseq_CUTTag_Find_meme_subparser.add_argument("-o","--output", required=True, type=str, help="Output directory")
    ChIPseq_CUTTag_Find_meme_subparser.add_argument("-t","--thread", default=10, type=int, help="Number of threads, default: 10")
    ChIPseq_CUTTag_Find_meme_subparser.add_argument("--nmotifs", default=3, type=int, help="Number of motifs to extract, default: 3")
    ChIPseq_CUTTag_Find_meme_subparser.add_argument("-db","--database", required=True, type=str, help="Database file")
    # homer heatmap
    ChIPseq_CUTTag_homer_heatmap_subparser = ChIPseq_CUTTag_motif_subparser.add_parser('homer_heatmap', help='Plot the homer finding result heatmap')
    ChIPseq_CUTTag_homer_heatmap_subparser.add_argument("-i","--input", required=True, type=str, help="The homer results information csv file")
    ChIPseq_CUTTag_homer_heatmap_subparser.add_argument("-o","--output", required=True, type=str, help="The output directory")
    ChIPseq_CUTTag_homer_heatmap_subparser.add_argument("--motif", default='no', type=str, help="The specified motifs txt file")
    ChIPseq_CUTTag_homer_heatmap_subparser.add_argument("--width", default=10, type=float, help="The plot width, default: 10")
    ChIPseq_CUTTag_homer_heatmap_subparser.add_argument("--height", default=8, type=float, help="The plot height, default: 8")
    # enrich heatmap
    ChIPseq_CUTTag_enrich_heatmap_subparser = ChIPseq_CUTTag_motif_subparser.add_parser('enrich_heatmap', help='Perform the motif enriched heatmap')
    ChIPseq_CUTTag_enrich_heatmap_subparser.add_argument("--motif_bed", required=True, type=str, help="The motif bed file information csv file")
    ChIPseq_CUTTag_enrich_heatmap_subparser.add_argument("--bed", required=True, type=str, help="The input bed file information csv file")
    ChIPseq_CUTTag_enrich_heatmap_subparser.add_argument("-o","--output", required=True, type=str, help="The output directory")
    ChIPseq_CUTTag_enrich_heatmap_subparser.add_argument("-t","--thread", default=12, type=int, help="The number of threads")
    ChIPseq_CUTTag_enrich_heatmap_subparser.add_argument("-s","--species", required=True, type=str, help="Species")

    # ChIPseq_CUTTag: 5. Accessbility
    ChIPseq_CUTTag_accessbility_parser = ChIPseq_CUTTag_subparsers.add_parser('accessbility', help='ChIPseq_CUTTag accessbility analysis')
    ChIPseq_CUTTag_accessbility_parser.add_argument("--bw", required=True, help="The input bw txt file")
    ChIPseq_CUTTag_accessbility_parser.add_argument("-o","--output", required=True, help="Output directory")
    ChIPseq_CUTTag_accessbility_parser.add_argument("--mode", required=True, choices=["region", "gene"], help="Analysis mode: region/gene")
    ChIPseq_CUTTag_accessbility_parser.add_argument("-i","--input", required=True, help="Input gene region or gene list txt file")
    ChIPseq_CUTTag_accessbility_parser.add_argument("--species", help="Species", choices=list(config['data'].keys()))
    ChIPseq_CUTTag_accessbility_parser.add_argument("--method", choices=["sample", "region_list"], default="sample", help="Boxplot plotting method：sample/region list")
    ChIPseq_CUTTag_accessbility_parser.add_argument("--center", required=True, help="Average profile alignment：region mode:peak/position，gene mode:TSS/TTS/whole")
    ChIPseq_CUTTag_accessbility_parser.add_argument("--promoter_size", type=int, default=5000, help="Gene promoter region size, default: TSS upstream and downstream 5kb")
    ChIPseq_CUTTag_accessbility_parser.add_argument("--region_list", default='', help="Region list file for box plot")

    # ChIPseq_CUTTag: 6. IGV snapshot
    ChIPseq_CUTTag_igv_parser = ChIPseq_CUTTag_subparsers.add_parser('igv', help='ChIPseq_CUTTag IGV snapshot')
    ChIPseq_CUTTag_igv_parser.add_argument("-i","--input", required=True, type=str, help="The input list file")
    ChIPseq_CUTTag_igv_parser.add_argument("-o","--output", required=True, type=str, help="Output directory")
    ChIPseq_CUTTag_igv_parser.add_argument("--genome", required=True, type=str, help="Select the genome id")
    ChIPseq_CUTTag_igv_parser.add_argument("--display_mode", default='collapse', help="Display mode", choices=["collapse","expand","squish"])
    ChIPseq_CUTTag_igv_parser.add_argument("--region", required=True, help="Specify the regions, separate by ,")
    ChIPseq_CUTTag_igv_parser.add_argument("--type", default='png', type=str, choices=['svg','png'], help="svg or png output")

    # ChIPseq_CUTTag: 7. Similarity
    ChIPseq_CUTTag_correlation_parser = ChIPseq_CUTTag_subparsers.add_parser('similarity', help='ChIPseq_CUTTag sample similarity analysis')
    ChIPseq_CUTTag_correlation_parser.add_argument("--bw1", required=True, help='Path to the first bw file')
    ChIPseq_CUTTag_correlation_parser.add_argument("--bw2", required=True, help='Path to the second bw file')
    ChIPseq_CUTTag_correlation_parser.add_argument("-o","--output", required=True, help='Output directory')
    ChIPseq_CUTTag_correlation_parser.add_argument("-t","--thread", type=int, default=30, help='Number of threads')
    ChIPseq_CUTTag_correlation_parser.add_argument("--region", help='Optional BED file for specific regions')

    return ChIPseq_CUTTag_parser
