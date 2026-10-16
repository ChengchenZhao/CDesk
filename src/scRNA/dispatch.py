"""
scRNA functions
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

# scRNA 1. Preprocess
def scRNA_preprocess(args):
    args.input = args.input.rstrip('/');args.output = args.output.rstrip('/')
    args.input = os.path.abspath(args.input);args.output = os.path.abspath(args.output)
    python_script = os.path.join(_root_dir,'src','scRNA','1.preprocess','preprocess.py')    
    if args.scRNA_preprocess_mode != 'smartseq':
        cmd = [_CDesk_env['python3_7'], python_script,args.scRNA_preprocess_mode,'-i',args.input,'-o',args.output,'-s',args.species,'-t',str(args.thread)]
        if args.scRNA_preprocess_mode == 'cellranger':
            cmd += ['--chemistry',args.chemistry]
        if args.scRNA_preprocess_mode == 'drseq':
            cmd += ['--umi',args.umi,'--barcode',args.barcode]
        subprocess.run(cmd,check=True,env=_CDesk_env)

    elif args.scRNA_preprocess_mode == 'smartseq':
        python_script = os.path.join(_root_dir,'1.RNAseq','1.Preprocess','data.py')
        cmd = [_CDesk_env['python3_7'], python_script,'-i',args.input,'-o',args.output]
        subprocess.run(cmd,check=True,env=_CDesk_env)
        bash_script = os.path.join(_root_dir,'1.RNAseq','1.Preprocess','BulkRNAseqData.sh')
        cmd = [bash_script,'-s', args.species, '-i', os.path.join(args.output,'data'), '-o', args.output, '-p', str(args.thread)]
        subprocess.run(cmd,check=True,env=_CDesk_env)

        if args.bw:
            bash_script = os.path.join(_root_dir,'1.RNAseq','1.Preprocess','BulkRNAseqBW.sh')
            cmd = [bash_script, '-i', f"{args.output}/Bam", '-o', f"{args.output}/BW", '-p', str(args.thread)]
            subprocess.run(cmd, check=True,env=_CDesk_env)
            
        if args.te:
            bash_script = os.path.join(_root_dir,'1.RNAseq','1.Preprocess','BulkRNAseqTE.sh')
            cmd = [bash_script, '-i', f"{args.output}/Bam", '-o', f"{args.output}/TE", '-p', str(args.thread),'-s',args.species]
            subprocess.run(cmd, check=True,env=_CDesk_env)

# scRNA 2. cluster
def scRNA_cluster(args):
    args.input = args.input.rstrip('/');args.output = args.output.rstrip('/')
    args.input = os.path.abspath(args.input);args.output = os.path.abspath(args.output)
    python_script = os.path.join(_root_dir,'src','scRNA','2.clustering','clustering.py')
    R_script = os.path.join(_root_dir,'src','scRNA','2.clustering','clustering.R')
    if args.mode == 'scanpy':
        cmd = [
        _CDesk_env['python3_7'], python_script,
        '-i', args.input,
        '-o', args.output,
        '-n', args.name,
        '--min_cells', str(args.min_cells),
        '--min_features', str(args.min_features),
        '--nFeature_RNA_min', str(args.nFeature_RNA_min),
        '--nFeature_RNA_max', str(args.nFeature_RNA_max),
        '--mt_percent', str(args.mt_percent),
        '--variable_features', str(args.variable_features),
        '--dim_prefer', str(args.dim_prefer),
        '--res', str(args.res),
        '--width',str(args.width),
        '--height',str(args.height)
        ]
        subprocess.run(cmd, check=True,env=_CDesk_env)

    elif args.mode == 'seurat':
        cmd = [
        'Rscript', R_script,
        args.input,
        args.output,
        args.name,
        str(args.min_cells),
        str(args.min_features),
        str(args.nFeature_RNA_min),
        str(args.nFeature_RNA_max),
        str(args.mt_percent),
        str(args.variable_features),
        str(args.dim_prefer),
        str(args.res),str(args.width),str(args.height)
        ]
        subprocess.run(cmd, check=True,env=_CDesk_env)

# scRNA: 3. Annotation
def scRNA_annotation(args):
    R_script = os.path.join(_root_dir,'src','scRNA','3.annotation','scRNAseqAnno.R')
    if args.scRNA_annotation_mode == 'marker':
        cmd = ['Rscript',R_script,'marker',args.input,args.marker,args.meta,args.output,str(args.expression_thre),
               str(args.percentage_thre),str(args.marker_thre),args.cluster,str(args.width),str(args.height)]
        subprocess.run(cmd, check=True,env=_CDesk_env)
    
    elif args.scRNA_annotation_mode == 'transfer':
        cmd = ['Rscript',R_script,'transfer',args.ref,args.query,args.meta,args.output,args.cluster,
               str(args.nfeatures),str(args.dims),str(args.width),str(args.height)]
        subprocess.run(cmd, check=True,env=_CDesk_env)

# scRNA 4. marker
def scRNA_marker(args):
    args.input = args.input.rstrip('/');args.output = args.output.rstrip('/')
    args.input = os.path.abspath(args.input);args.output = os.path.abspath(args.output)
    if args.scRNA_marker_mode == 'all':
        R_script = os.path.join(_root_dir,'src','scRNA','4.marker','findAllMarkers.R')
        cmd = [
            "Rscript", R_script,
            args.input, args.meta,str(args.fc), str(args.p),str(args.m), args.output,str(args.width),str(args.height)
        ]
        subprocess.run(cmd, check=True,env=_CDesk_env)

    elif args.scRNA_marker_mode == '2group':
        R_script = os.path.join(_root_dir,'src','scRNA','4.marker','marker_in_2_group.R')
        cmd = [
            "Rscript", R_script,
            '-i',args.input,'-g',args.meta ,'-t',args.type1,'-e',args.type2,
            '-f',str(args.fc), '--adjPvalFilter',str(args.p),'-m',str(args.m),'-o',args.output,
            '--width',str(args.width),'--height',str(args.height)
        ]
        subprocess.run(cmd, check=True,env=_CDesk_env)

# scRNA: 5. trajectory
def scRNA_trajectory(args):
    if args.scRNA_trajectory_mode == 'monocle':
        args.input = args.input.rstrip('/');args.output = args.output.rstrip('/')
        args.input = os.path.abspath(args.input);args.output = os.path.abspath(args.output)
        R_script = os.path.join(_root_dir,'src','scRNA','5.trajectory','monocle3.R')
        cmd = ['Rscript',R_script,args.input,args.output,args.meta,args.start,str(args.width),str(args.height)]
        subprocess.run(cmd, check=True,env=_CDesk_env)

    elif args.scRNA_trajectory_mode == 'diffusion':
        args.input = args.input.rstrip('/');args.output = args.output.rstrip('/')
        args.input = os.path.abspath(args.input);args.output = os.path.abspath(args.output)
        R_script = os.path.join(_root_dir,'src','scRNA','5.trajectory','diffusionmap.R')
        cmd = ['Rscript',R_script,args.input,args.output,args.meta,args.start,str(args.width),str(args.height)]
        subprocess.run(cmd, check=True,env=_CDesk_env)

    elif args.scRNA_trajectory_mode == 'velocity':
        args.cellranger = args.cellranger.rstrip('/');args.output = args.output.rstrip('/');args.rds = args.rds.rstrip('/');args.genes = args.genes.rstrip('/')
        args.cellranger = os.path.abspath(args.cellranger);args.output = os.path.abspath(args.output);args.rds = os.path.abspath(args.rds);args.genes = os.path.abspath(args.genes) 
        python_script = os.path.join(_root_dir,'src','scRNA','5.trajectory','velocity.py')
        cmd = [_CDesk_env['python3_7'],python_script,'--cellranger_output',args.cellranger,
               '--rds_file',args.rds,'--gtf',args.gtf,'--genes',args.genes,'-o',args.output,
               '--meta',args.meta,'--width',str(args.width),'--height',str(args.height),'-t',str(args.thread)]
        subprocess.run(cmd, check=True,env=_CDesk_env)
        figures_dir = os.path.join(os.getcwd(),'figures')
        if os.path.exists(figures_dir) and os.path.isdir(figures_dir):
            shutil.rmtree(figures_dir)

    elif args.scRNA_trajectory_mode == 'cstreet':
        args.expression = args.expression.rstrip('/');args.output = args.output.rstrip('/');args.state = args.state.rstrip('/')
        args.expression = os.path.abspath(args.expression);args.output = os.path.abspath(args.output);args.state = os.path.abspath(args.state)
        python_script = os.path.join(_root_dir,'src','scRNA','5.trajectory','cstreet.py')
        cmd = [_CDesk_env['python3_7'],python_script,'-i',args.expression,'-s',args.state,'-n',args.name,'-o',args.output]
        subprocess.run(cmd, check=True,env=_CDesk_env)

# scRNA: 6. similarity
def scRNA_similarity(args):
    args.output = args.output.rstrip('/');args.sample1 = args.sample1.rstrip('/');args.sample2 = args.sample2.rstrip('/');args.type1 = args.type1.rstrip('/');args.type2 = args.type2.rstrip('/')
    args.output = os.path.abspath(args.output);args.sample1 = os.path.abspath(args.sample1);args.sample2 = os.path.abspath(args.sample2);args.type1 = os.path.abspath(args.type1);args.type2 = os.path.abspath(args.type2)
    R_script = os.path.join(_root_dir,'src','scRNA','6.similarity','correlation.R')
    cmd = ['Rscript',R_script,args.sample1,args.sample2,args.meta1,args.meta2,args.type1,args.type2,args.output,str(args.nfeatures_used),str(args.integrate_features),str(args.width),str(args.height)]
    subprocess.run(cmd, check=True,env=_CDesk_env)

# scRNA: 7. interaction
def scRNA_interaction(args):
    args.input = args.input.rstrip('/');args.output = args.output.rstrip('/')
    args.input = os.path.abspath(args.input);args.output = os.path.abspath(args.output)
    if args.scRNA_interaction_mode == 'cellchat':
        R_script = os.path.join(_root_dir,'src','scRNA','7.interaction','cellchat.R')
        cmd = ['Rscript',R_script,'--seurat_file',args.input,"--output_directory",args.output,"--threads",str(args.thread),"--species",args.species,"--group_by",args.group,"--width",str(args.width),"--height",str(args.height)]
        subprocess.run(cmd, check=True,env=_CDesk_env)

    elif args.scRNA_interaction_mode == 'pyscenic':
        args.tf = args.tf.rstrip('/');args.db = args.db.rstrip('/')
        args.tf = os.path.abspath(args.tf);args.db = os.path.abspath(args.db)
        python_script = os.path.join(_root_dir,'src','scRNA','7.interaction','GRN.py')
        cmd = [_CDesk_env['python3_7'],python_script,'-i',args.input,'-o',args.output,"--tf-list",args.tf,"--db-folder",args.db,"--species",args.species,"--num-workers",str(args.thread),"--auc-threshold",str(args.auc),"--nes-threshold",str(args.nes),'--normalize',str(args.normalize),'--scale',str(args.scale),'--log',str(args.log),'--top_regulons',str(args.top_regulons),'--top_targets',str(args.top_targets),'--plot_tf',str(args.plot_tf)]
        subprocess.run(cmd, check=True,env=_CDesk_env)

# scRNA: 8. integrate
def scRNA_integrate(args):
    if args.scRNA_integrate_mode == 'seurat':
        args.input = args.input.rstrip('/');args.output = args.output.rstrip('/')
        args.input = os.path.abspath(args.input);args.output = os.path.abspath(args.output)
        R_script = os.path.join(_root_dir,'src','scRNA','8.integrate','integrate.R')
        cmd = ['Rscript',R_script,args.input,args.mode,args.output,str(args.min_cells),str(args.min_features),str(args.nfeatures_used),
               str(args.integrate_features),str(args.pc_num),str(args.pc_reduce),str(args.max_iter_harmony),str(args.k_anchor),str(args.width),str(args.height)]
        subprocess.run(cmd, check=True,env=_CDesk_env)

    elif args.scRNA_integrate_mode == 'bulk':
        args.meta = args.meta.rstrip('/');args.output = args.output.rstrip('/');args.bulk = args.bulk.rstrip('/');args.scRNA = args.scRNA.rstrip('/')
        args.meta = os.path.abspath(args.meta);args.output = os.path.abspath(args.output);args.bulk = os.path.abspath(args.bulk);args.scRNA = os.path.abspath(args.scRNA)
        R_script = os.path.join(_root_dir,'src','scRNA','8.integrate','bulk+scRNA.R')
        cmd = ['Rscript',R_script,args.bulk,args.scRNA,args.meta,args.output,str(args.nfeatures_used),
               str(args.integrate_features),str(args.pc_num),str(args.pc_reduce),str(args.width),str(args.height)]
        subprocess.run(cmd, check=True,env=_CDesk_env)

def handle(args):
    """scRNA command dispatcher"""
    cmd = args.scRNA_command
    if cmd == 'preprocess':
        scRNA_preprocess(args)
    elif cmd == 'cluster':
        scRNA_cluster(args)
    elif cmd == 'annotation':
        scRNA_annotation(args)
    elif cmd == 'marker':
        scRNA_marker(args)
    elif cmd == 'trajectory':
        scRNA_trajectory(args)
    elif cmd == 'similarity':
        scRNA_similarity(args)
    elif cmd == 'interaction':
        scRNA_interaction(args)
    elif cmd == 'integrate':
        scRNA_integrate(args)
