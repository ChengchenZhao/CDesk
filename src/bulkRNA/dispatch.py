"""
bulkRNA functions
"""

import os
import sys
import subprocess
import shutil

_root_dir = None
_CDesk_env = None
_config = None

def init(root_dir, CDesk_env, config=None):
    global _root_dir, _CDesk_env, _config
    _root_dir = root_dir
    _CDesk_env = CDesk_env
    _config = config

# bulkRNA 1. preprocess
def bulkRNA_preprocess(args): 
    args.input = args.input.rstrip('/');args.output = args.output.rstrip('/')
    args.input = os.path.abspath(args.input);args.output = os.path.abspath(args.output)
    python_script = os.path.join(_root_dir,'src','bulkRNA','1.Preprocess','data.py')
    cmd = [_CDesk_env['python3_7'], python_script,'-i',args.input,'-o',args.output]
    subprocess.run(cmd,check=True,env=_CDesk_env)

    bash_script = os.path.join(_root_dir,'src','bulkRNA','1.Preprocess','BulkRNAseqData.sh')
    cmd = [bash_script,'-s', args.species, '-i', os.path.join(args.output,'data'), '-o', args.output, '-p', str(args.thread)]
    subprocess.run(cmd,check=True,env=_CDesk_env)

    if args.bw:
        bash_script = os.path.join(_root_dir,'src','bulkRNA','1.Preprocess','BulkRNAseqBW.sh')
        cmd = [bash_script, '-i', f"{args.output}/Bam", '-o', f"{args.output}/BW", '-p', str(args.thread)]
        subprocess.run(cmd, check=True,env=_CDesk_env)

    if args.te:
        bash_script = os.path.join(_root_dir,'src','bulkRNA','1.Preprocess','BulkRNAseqTE.sh')
        cmd = [bash_script, '-i', f"{args.output}/Bam", '-o', f"{args.output}/TE", '-p', str(args.thread),'-s',args.species]
        subprocess.run(cmd, check=True,env=_CDesk_env)

# bulkRNA 2. cor
def bulkRNA_cor(args):
    R_script = os.path.join(_root_dir,'src','bulkRNA','2.correlation','cor.R')
    args.input = args.input.rstrip('/');args.output = args.output.rstrip('/');args.group = args.group.rstrip('/')
    args.input = os.path.abspath(args.input);args.output = os.path.abspath(args.output);args.group = os.path.abspath(args.group)
    cmd = ['Rscript', R_script,args.input,args.output,args.group,str(args.batch),str(args.width),str(args.height)]
    subprocess.run(cmd, check=True,env=_CDesk_env)

# bulkRNA 3. DEG
def bulkRNA_DEG(args):
    python_script = os.path.join(_root_dir,'src','bulkRNA','3.DEGD','gfold.py')
    args.input = args.input.rstrip('/');args.output = args.output.rstrip('/')
    args.input = os.path.abspath(args.input);args.output = os.path.abspath(args.output)
    if args.DEG_mode != 'gfold':
        args.group = args.group.rstrip('/');args.group = os.path.abspath(args.group)
        R_script = os.path.join(_root_dir,'src','bulkRNA','3.DEGD','degd.R')
        cmd = [
            'Rscript', R_script,args.DEG_mode,
            args.input, args.output,args.group,str(args.p),args.gene,str(args.fc),
            str(args.pval),str(args.top_num),str(args.width),str(args.height)
        ]
        subprocess.run(cmd, check=True,env=_CDesk_env)

    elif args.DEG_mode == 'gfold':
        cmd = [
                _CDesk_env['python3_7'], python_script,
                "--species",args.species,"--output_dir",args.output,
                "--meta_file",args.input,"--top_genes",str(args.top_num),
                "--gene_of_interest",args.gene,"--fc_threshold",str(args.gfold),
                "--thread",str(args.thread),"--width",str(args.width),"--height",str(args.height)
        ]

        if args.p:
            cmd = cmd + ['--plot']

        subprocess.run(cmd, check=True,env=_CDesk_env)

# bulkRNA 4. enrich analysis
def bulkRNA_enrich(args):
    if args.enrich_mode == 'analyze':
        args.input = args.input.rstrip('/');args.output = args.output.rstrip('/')
        args.input = os.path.abspath(args.input);args.output = os.path.abspath(args.output)
        R_script = os.path.join(_root_dir,'src','bulkRNA','4.enrich','enrichment_analysis.R')
        if args.type == 'multi':
            R_script = os.path.join(_root_dir,'src','bulkRNA','4.enrich','enrichment_multi.R')
        if args.species == 'no' and args.custom == 'no':
            print('At least provide species or custom file!')
            sys.exit(1)
        if args.custom != 'no':
            args.custom = args.custom.rstrip('/');args.custom = os.path.abspath(args.custom)
        cmd = ["Rscript", R_script,args.input,args.custom,args.species,args.output,str(args.width),str(args.height)]
        subprocess.run(cmd, check=True,env=_CDesk_env)

    elif args.enrich_mode == 'plot':
        R_script = os.path.join(_root_dir,'src','bulkRNA','4.enrich','enrichment_plot.R')
        if args.type == 'multi':
            R_script = os.path.join(_root_dir,'src','bulkRNA','4.enrich','enrichment_multi_plot.R')
        args.input = args.input.rstrip('/');args.output = args.output.rstrip('/')
        args.input = os.path.abspath(args.input);args.output = os.path.abspath(args.output)
        cmd = ["Rscript", R_script, args.input,args.output,str(args.width),str(args.height)]
        subprocess.run(cmd, check=True,env=_CDesk_env)

    elif args.enrich_mode == 'GSEA': 
        R_script = os.path.join(_root_dir,'src','bulkRNA','4.enrich','gsea.R')
        args.input = args.input.rstrip('/');args.output = args.output.rstrip('/')
        args.input = os.path.abspath(args.input);args.output = os.path.abspath(args.output)
        cmd = ['Rscript',R_script,args.input,args.path,args.species,args.cols,args.output,str(args.width),str(args.height)]
        subprocess.run(cmd, check=True,env=_CDesk_env)

# bulkRNA 5. Similarity Analysis
def bulkRNA_Similarity(args):
    R_script = os.path.join(_root_dir,'src','bulkRNA','5.SimilarityAnalysis','Similarity.R')
    args.output = args.output.rstrip('/');args.input = args.input.rstrip('/');args.group = args.group.rstrip('/')
    args.output = os.path.abspath(args.output);args.input = os.path.abspath(args.input);args.group = os.path.abspath(args.group)
    cmd = ['Rscript',R_script,args.input,args.group,args.output,args.gene,args.batch,str(args.width),str(args.height)]
    subprocess.run(cmd, check=True,env=_CDesk_env)

# bulkRNA: 6. Clustering
def bulkRNA_Clustering(args):
    if args.cluster_mode == 'WGCNA':
        R_script = os.path.join(_root_dir,'src','bulkRNA','6.Clustering','WGCNA.R')
        args.output = args.output.rstrip('/');args.input = args.input.rstrip('/');args.pheno = args.pheno.rstrip('/')
        args.output = os.path.abspath(args.output);args.input = os.path.abspath(args.input);args.pheno = os.path.abspath(args.pheno)
        cmd = ['Rscript', R_script,args.input,args.pheno,args.trait,args.output]
        subprocess.run(cmd, check=True,env=_CDesk_env)
    elif args.cluster_mode == 'kmean':
        R_script = os.path.join(_root_dir,'src','bulkRNA','6.Clustering','Kmeans.R')
        args.output = args.output.rstrip('/'); args.input = args.input.rstrip('/');args.group = args.group.rstrip('/')
        args.output = os.path.abspath(args.output);args.input = os.path.abspath(args.input);args.group = os.path.abspath(args.group)
        if args.gene != 'no':
            args.gene = args.gene.rstrip('/');args.gene = os.path.abspath(args.gene)
        cmd = ['Rscript', R_script,args.input,args.group,args.output,str(args.cluster),args.method,
               str(args.height),str(args.width),args.gene]
        subprocess.run(cmd, check=True,env=_CDesk_env)

# bulkRNA 7. AlternativeSplicingDetection
def bulkRNA_AlternativeSpliceDetection(args):
    args.output = args.output.rstrip('/');args.b1 = args.b1.rstrip('/')
    args.output = os.path.abspath(args.output);args.b1 = os.path.abspath(args.b1)
    if args.b2 != 'No':
        args.b2 = args.b2.rstrip('/');args.b2 = os.path.abspath(args.b2)
    if args.Splice_mode == 'detect':
        python_script = os.path.join(_root_dir,'src','bulkRNA','7.AlternativeSplicingDetection','alternative_splicing.py')
        cmd = [_CDesk_env['python3_6'], python_script,'--b1',args.b1,'--b2',args.b2,'-s',args.species,'-o',args.output,
                "--type",args.type,"-t",str(args.thread),'--length',str(args.length),
                '--variable_read_length',str(args.variable_read_length),'--allow_clipping',str(args.allow_clipping)]
        subprocess.run(cmd, check=True,env=_CDesk_env)
        
    elif args.Splice_mode == 'draw':
        python_script = os.path.join(_root_dir,'src','bulkRNA','7.AlternativeSplicingDetection','draw_sashimiplot.py')
        cmd = [_CDesk_env['python3_7'], python_script,'--b1',args.b1,'--b2',args.b2,'-s',args.species,'-o',args.output,
                "--region",args.region,"--group",args.group,'--l1',args.l1,'--l2',args.l2,
                '--exon_s',str(args.exon_s),'--intron_s',str(args.intron_s)]
        subprocess.run(cmd,check=True,env=_CDesk_env)

def handle(args):
    cmd = args.bulkRNA_command
    if cmd == 'preprocess':
        bulkRNA_preprocess(args)
    elif cmd == 'correlation':
        bulkRNA_cor(args)
    elif cmd == 'DEG':
        bulkRNA_DEG(args)
    elif cmd == 'enrich':
        bulkRNA_enrich(args)
    elif cmd == 'similarity':
        bulkRNA_Similarity(args)
    elif cmd == 'cluster':
        bulkRNA_Clustering(args)
    elif cmd == 'splice':
        bulkRNA_AlternativeSpliceDetection(args)