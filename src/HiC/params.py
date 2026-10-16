"""
HiC params
"""

def register_subparser(subparsers, config):
    HiC_parser = subparsers.add_parser('HiC', help='HiC pipeline')
    HiC_subparsers = HiC_parser.add_subparsers(
        help='HiC pipeline', title='HiC pipeline', dest='HiC_command'
    )

    # HiC: 1. Preprocess
    HiC_preprocess_parser = HiC_subparsers.add_parser('preprocess', help='HiC preprocess pipeline')
    HiC_preprocess_parser.add_argument("-i","--input", required=True, type=str, help="Input fastq file information file")
    HiC_preprocess_parser.add_argument("-o","--output", required=True, type=str, help="Output directory")
    HiC_preprocess_parser.add_argument('-s','--species', required=True, type=str, help='Specify the species', choices=list(config['data'].keys()))
    HiC_preprocess_parser.add_argument('-t','--thread', type=int, default=10, help="Number of thread, default:10")
    HiC_preprocess_parser.add_argument('-f','--filter', action='store_true', help="Do stricter filtering, use parameters refering to Ruan's Lab")
    HiC_preprocess_parser.add_argument('-r','--resolution', type=str, default='10000,100000,1000000', help="Specify the resolutions, separate by comma, default: 10000,100000,1000000")

    # HiC: 2. Correlation
    HiC_correlation_parser = HiC_subparsers.add_parser('correlation', help='HiC sample correlation analysis')
    HiC_correlation_parser.add_argument('-i',"--input", type=str, required=True, help="The HiC input sample information file")
    HiC_correlation_parser.add_argument('-o','--output', type=str, required=True, help="Specify the output directory")
    HiC_correlation_parser.add_argument('--chr', type=str, default='No', help="Specify the chromosome")
    HiC_correlation_parser.add_argument('--inter_chrom', action='store_true', help="Whether to include inter chrom signals")
    HiC_correlation_parser.add_argument('--width', type=float, default=8, help="The width of the plot, default:8")
    HiC_correlation_parser.add_argument('--height', type=float, default=6, help="The height of the plot, default:6")

    # HiC: 3. Matrix
    HiC_matrix_parser = HiC_subparsers.add_parser('matrix', help='HiC matrix balancing and format tranfer')
    HiC_matrix_parser.add_argument('--inputFormat', type=str, required=True, help="Input format", choices=['pairs','cool','mcool','juicerhic','fanchic','h5'])
    HiC_matrix_parser.add_argument('--outputFormat', type=str, required=True, help="Output format", choices=['cool','mcool','juicerhic','fanchic','h5'])
    HiC_matrix_parser.add_argument('-i','--input', type=str, required=True, help="Input file")
    HiC_matrix_parser.add_argument('-o','--output', type=str, required=True, help="Output directory")
    HiC_matrix_parser.add_argument('-t','--thread', type=int, default=10, help="Number of threads, default: 10")
    HiC_matrix_parser.add_argument('-r','--resolution', type=str, default='', help="Resolutions")
    HiC_matrix_parser.add_argument('-n','--norm', type=str, default='None', help="Norm method", choices=['KR','ICE','VC','VC-SQRT','None'])
    HiC_matrix_parser.add_argument('--oe', type=str, default='', help="The resolution of O/E matrix to get")
    HiC_matrix_parser.add_argument('--sparse', type=str, default='', help="The resolution of sparse matrix to get")
    HiC_matrix_parser.add_argument('-s','--species', default='', type=str, help='Specify the species', choices=list(config['data'].keys()))

    # HiC: 4. TAD
    HiC_tad_parser = HiC_subparsers.add_parser('TAD', help='HiC TAD analysis')
    HiC_tad_parser.add_argument('-i', '--input', type=str, required=True, help="The information file")
    HiC_tad_parser.add_argument('-o', '--output', type=str, required=True, help="The output directory")
    HiC_tad_parser.add_argument('--insulation_window', type=int, default=1000000, help="Window sizes to calculate insulation score in base pairs, default: 1000000")
    HiC_tad_parser.add_argument('--directionality_window', type=int, default=1000000, help="Window sizes to calculate directionality index in base pairs, default: 1000000")
    HiC_tad_parser.add_argument('--expand', type=int, default=1000000, help="Window sizes of boundary expansion in base pairs to calulate insulation score and DI near/on boundary, default: 1000000")

    # HiC: 5. Compartment
    HiC_compartment_parser = HiC_subparsers.add_parser('compartment', help='HiC compartment analysis')
    HiC_compartment_parser.add_argument('-i', '--input', type=str, required=True, help="The information file")
    HiC_compartment_parser.add_argument('-r', '--resolution', type=int, required=True, help="Resolution to specify")
    HiC_compartment_parser.add_argument('-o', '--output', type=str, required=True, help="The output directory")
    HiC_compartment_parser.add_argument('-t', '--thread', type=int, default=20, help="The number of threads, default: 20")
    HiC_compartment_parser.add_argument('--GC', type=str, default='no', help="The GC file to provide, default: no")
    HiC_compartment_parser.add_argument('-s', '--species', type=str, default='no', help="The species, need to specify if no GC file provided, default: no", choices=['no'] + list(config['data'].keys()))
    HiC_compartment_parser.add_argument('--width', type=float, default=10, help="Default: 10")
    HiC_compartment_parser.add_argument('--height', type=float, default=8, help="Default: 8")

    # HiC: 6. Loop
    HiC_loop_parser = HiC_subparsers.add_parser('loop', help='HiC loop calling')
    HiC_loop_parser.add_argument("-i","--input", required=True, type=str, help="FanC supported hic input file")
    HiC_loop_parser.add_argument("-o","--output", required=True, type=str, help="Output directory")
    HiC_loop_parser.add_argument('-t','--thread', type=int, default=10, help="Number of thread, default:10")

    # HiC: 7. Reconstruction
    HiC_reconstruction_parser = HiC_subparsers.add_parser('reconstruct', help='HiC 3D reconstuction')
    # FLAMINGO
    HiC_reconstruction_parser.add_argument('-i','--input', required=True, help='The input hic file')
    HiC_reconstruction_parser.add_argument("-o","--output", required=True, type=str, help="Output directory")
    HiC_reconstruction_parser.add_argument('--chr', required=True, help='The chromosome')
    HiC_reconstruction_parser.add_argument('--domain_res', required=True, help='The domain resolution')
    HiC_reconstruction_parser.add_argument('--frag_res', required=True, help='The fragment resolution')
    HiC_reconstruction_parser.add_argument('--norm', default='NONE', help='The normalization method, default: NONE')
    HiC_reconstruction_parser.add_argument('-t','--thread', type=int, default=10, help="Number of thread, default: 10")

    # HiC: 8. pss
    HiC_pss_parser = HiC_subparsers.add_parser('pss', help='HiC distance contact analysis and pss plot')
    HiC_pss_parser.add_argument("-i","--input", required=True, type=str, help="Input sample information file")
    HiC_pss_parser.add_argument("-o", "--output", required=True, type=str, help="Output directory")
    HiC_pss_parser.add_argument('--bin', type=str, default='20', help="Bins divided between each order of magnitude, default: 20")
    HiC_pss_parser.add_argument('--curve_mean_min', type=str, default='1e-6', help="The lower bound of the mean of the curve, default: 1e-6")
    HiC_pss_parser.add_argument('--curve_bin_min', type=str, default='0', help="The lower bound of bin index of the curve, deafault: 0")
    HiC_pss_parser.add_argument('--curve_x_min', type=str, default='4', help="The minimal X-axis truncation value of the curve, default: 4")
    HiC_pss_parser.add_argument('--curve_x_max', type=str, default='8', help="The maximal X-axis truncation value of the curve, default: 8")
    HiC_pss_parser.add_argument('--expcurve_y_bias', type=str, default='-4', help="The Y-axis offset when plotting the difference from the expectation curve, default: -4")
    HiC_pss_parser.add_argument('--expcurve_x_min', type=str, default='5', help="The minimal X-axis drop value when plotting the difference from the expectation curve, default: 5")
    HiC_pss_parser.add_argument('--expcurve_x_max', type=str, default='8', help="The maximal X-axis drop value when plotting the difference from the expectation curve, default: 8")
    HiC_pss_parser.add_argument('--heatmap_x_min', type=str, default='4.8', help="The minimal X-axis truncation value of the heatmap, default: 4.8")
    HiC_pss_parser.add_argument('--heatmap_x_max', type=str, default='7.5', help="The maximal X-axis truncation value of the heatmap, default: 7.5")
    HiC_pss_parser.add_argument('--width', type=str, default='10', help="The width of the plot, default: 10")
    HiC_pss_parser.add_argument('--height', type=str, default='6', help="The height of the plot, default: 6")

    # HiC: 9. Compare
    HiC_compare_parser = HiC_subparsers.add_parser('compare', help='HiC contact compare analysis')
    HiC_compare_parser.add_argument('--hic1', type=str, required=True, help="The first hic file")
    HiC_compare_parser.add_argument('--hic2', type=str, required=True, help="The second hic file")
    HiC_compare_parser.add_argument('-o','--output', type=str, required=True, help="Specify the output directory")
    HiC_compare_parser.add_argument('--chr', type=str, required=True, help="Specify the chromosome")
    HiC_compare_parser.add_argument('--A_min', type=float, default=None, help="The average expression to filter out interactions, default: None")

    return HiC_parser
