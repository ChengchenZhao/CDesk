import os
import argparse
import sys
import glob
import subprocess

def cellranger(args):
    fastq_1_files = glob.glob(os.path.join(args.i, '*_1.fastq.gz'))
    fastq_2_files = glob.glob(os.path.join(args.i, '*_2.fastq.gz'))
    # 检查文件是否存在
    if not fastq_1_files:
        print("No _1.fastq.gz files found.")
        sys.exit(1)
    if not fastq_2_files:
        print("No _2.fastq.gz files found.")
        sys.exit(1)
    
    if not os.path.exists(args.o):
        os.makedirs(args.o)
    
    # cellranger
    transcriptome = {    
        "pig": "/mnt/liudong/data/Genome/susScr11Self",
        "test": "/mnt/linzejie/scRNA/count/refdata-gex-GRCh38-2020-A"
    }

    #if not os.path.exists(f'{args.o}/cellranger'):
    os.makedirs(f'{args.o}/cellranger')
    #if not os.path.exists(f'{args.o}/cellranger/input'):
    os.makedirs(f'{args.o}/cellranger/input')

    prefixes = list()
    # 软链接，修改文件名
    for file_name in os.listdir(args.i):
        if file_name.endswith('_1.fastq.gz'):
            base_name = file_name.replace('_1.fastq.gz', '')
            prefixes.append(base_name)
            symlink_name = f"{base_name}_S1_L001_R1_001.fastq.gz"
            source_file = os.path.join(args.i, file_name)
            target_file = os.path.join(f'{args.o}/cellranger/input', symlink_name)
            os.symlink(source_file, target_file)

        elif file_name.endswith('_2.fastq.gz'):
            base_name = file_name.replace('_2.fastq.gz', '')
            symlink_name = f"{base_name}_S1_L001_R2_001.fastq.gz"

            source_file = os.path.join(args.i, file_name)
            target_file = os.path.join(f'{args.o}/cellranger/input', symlink_name)

            os.symlink(source_file, target_file)

    # 运行cellranger
    for sample in prefixes:
        cmd = (f"cellranger count --id={sample} --fastqs={args.o}/cellranger/input --sample={sample} --transcriptome={transcriptome[args.s]} --localcores {args.t}")
        os.system(cmd)
        os.system(f'mv {sample} {args.o}/cellranger')

def drseq(args):
    fastq_1_files = glob.glob(os.path.join(args.i, '*_1.fastq'))
    fastq_2_files = glob.glob(os.path.join(args.i, '*_2.fastq'))
    # 检查文件是否存在
    if not fastq_1_files:
        print("No _1.fastq files found.")
        sys.exit(1)
    if not fastq_2_files:
        print("No _2.fastq files found.")
        sys.exit(1)
    
    if not os.path.exists(args.o):
        os.makedirs(args.o)

    prefixes = list()
    for file_name in os.listdir(args.i):
        if file_name.endswith('_1.fastq'):
            base_name = file_name.replace('_1.fastq', '')
            prefixes.append(base_name)

    # DrSeq
    refgene = {    
        "pig": "/mnt/liudong/data/Genome/susScr11/susScr11.DrSeq.refgenes.txt"
    }
    mapindex = {
        "pig":"/mnt/linzejie/temp/test/susScr11_bowtie2/susScr11"
    }

    #if not os.path.exists(f'{args.o}/DrSeq'):
    os.makedirs(f'{args.o}/DrSeq')
    # 运行DrSeq
    for sample in prefixes:
        cmd = (f"DrSeq simple -b {args.i}/{sample}_1.fastq -r {args.i}/{sample}_2.fastq -n {sample} -g {refgene[args.s]} --maptool bowtie2 --mapindex {mapindex[args.s]} --thread {args.t} -f")
        os.system(cmd)
        os.system(f'mv {sample} {args.o}/DrSeq')
        os.system(f'mv {sample}.conf {args.o}/DrSeq')

def preprocess_smartseq(args):
    bash_script = '../../1.RNAseq/1.PreProcessPipeline/BulkRNAseqData.sh'
    cmd = (f"{bash_script} -s {args.s} -i {args.i} -o {args.o} -p {args.t} -l {args.l}")
    os.system(cmd)
    if args.bw:
        bash_script = '../../1.RNAseq/1.PreProcessPipeline/BulkRNAseqBW.sh'
        cmd = (f"{bash_script} -i {args.i} -o {args.o}/BW -p {args.t}")
        os.system(cmd)
    if args.te:
        bash_script = '../../1.RNAseq/1.PreProcessPipeline/BulkRNAseqTE.sh'
        cmd = (f"{bash_script} -i {args.i} -o {args.o}/TE -p {args.t}")
        os.system(cmd)

def preprocess_singleron(args):
    command = f"conda run -n celescope python -c 'import sys; print(sys.executable)'"
    subprocess.run(command, shell=True, check=True)
    print("激活celescope Conda 环境")
    fastq_1_files = glob.glob(os.path.join(args.i, '*_1.fastq.gz'))
    fastq_2_files = glob.glob(os.path.join(args.i, '*_2.fastq.gz'))
    # 检查文件是否存在
    if not fastq_1_files:
        print("No _1.fastq.gz files found.")
        sys.exit(1)
    if not fastq_2_files:
        print("No _2.fastq.gz files found.")
        sys.exit(1)
    
    if not os.path.exists(args.o):
        os.makedirs(args.o)

    # singleron
    mapindex = {    
        "galGal": "/mnt/linzejie/temp/galGal6_celescope"
    }
    #if not os.path.exists(f'{args.o}/singleron'):
    os.makedirs(f'{args.o}/singleron')

    prefixes = list()
    for file_name in os.listdir(args.i):
        if file_name.endswith('_1.fastq.gz'):
            base_name = file_name.replace('_1.fastq.gz', '')
            prefixes.append(base_name)

    filename = f'{args.o}/singleron/rna.mapfile'
    with open(filename, 'w') as f:
        lines = [f'{x}\t{args.i}\t{x}' for x in prefixes]
        f.write('\n'.join(lines))

    #cmd = (f"multi_rna --mapfile {args.o}/singleron/rna.mapfile --genomeDir {mapindex[args.s]} --thread {args.t} --mod shell --outdir {args.o}/singleron")
    #os.system(cmd)
    cmd = ['conda', 'run', '-n', 'celescope', 'multi_rna', '--mapfile', f"{args.o}/singleron/rna.mapfile",
       '--genomeDir', mapindex[args.s],
       '--thread', str(args.t),
       '--mod', 'shell',
       '--outdir', f"{args.o}/singleron"]
    subprocess.run(cmd, check=True)
    os.system(f'mv shell {args.o}/singleron')

    for file_name in os.listdir(f'{args.o}/singleron/shell'):
        base_name = file_name.replace('_.sh', '')
        #cmd = (f'bash {args.o}/singleron/shell/{file_name}')
        #cmd = ['bash', f"{args.o}/singleron/shell/{file_name}"]
        cmd = ['conda', 'run', '-n', 'celescope', 'bash', f"{args.o}/singleron/shell/{file_name}"]
        subprocess.run(cmd, check=True)

def main():
    # CDesk ArgumentParser Object
    parser = argparse.ArgumentParser(prog='scRNA', description="scRNA pipeline")
    # 创建二级子命令
    subparsers = parser.add_subparsers(dest='type', required=True, help="Choose an omic to run the script")

    # 10x
    # 创建cellranger子命令的子解析器
    cellranger_parser = subparsers.add_parser('cellranger', help='cellranger 10x preprocess')
    cellranger_parser.add_argument('-i', required=True, help='Specify the input fastq file directory, fastqs file format: X_1.fastq.gz, X_2.fastq.gz')
    cellranger_parser.add_argument('-s', required=True, help='Specify the species')
    cellranger_parser.add_argument('-o', required=True, help='Specify the absolute directory of output file, the ending does not require a /')
    cellranger_parser.add_argument('-t', help='Specify number of threads (default is 8)',default=8,type=int)

    # 创建DrSeq子命令的子解析器
    drseq_parser = subparsers.add_parser('drseq', help='DrSeq 10x preprocess')
    drseq_parser.add_argument('-i', required=True, help='Specify the input fastq file directory, fastqs file format: X_1.fastq, X_2.fastq')
    drseq_parser.add_argument('-s', required=True, help='Specify the species')
    drseq_parser.add_argument('-o', required=True, help='Specify the absolute directory of output file, the ending does not require a /')
    drseq_parser.add_argument('-t', help='Specify number of threads (default is 8)',default=8,type=int)

    # singleron
    singleron_parser = subparsers.add_parser('singleron',help='singleron pipeline')
    singleron_parser.add_argument('-i', required=True, help='Specify the input fastq file directory, fastqs file format: X_1.fastq.gz, X_2.fastq.gz')
    singleron_parser.add_argument('-o', required=True, help='Specify the absolute directory of output file, the ending does not require a /')
    singleron_parser.add_argument('-s', required=True, help='Specify the species',
                                      choices=["galGal","pig"])
    singleron_parser.add_argument('-t', help='Specify number of threads (default is 8)',default=8,type=int)

    # 解析命令行参数
    args = parser.parse_args()

    # 判断并执行 bulkRNA 相关子命令
    if args.type == 'cellranger':
        cellranger(args)
    elif args.type == 'drseq':
        drseq(args)
    elif args.type == 'singleron':
        preprocess_singleron(args)

if __name__ == "__main__":
    main()
