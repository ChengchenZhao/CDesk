"""
ATAC functions
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

# ATAC: 1. Preprocess
def ATAC_preprocess(args):

    args.input = args.input.rstrip('/');args.output = args.output.rstrip('/')

    args.input = os.path.abspath(args.input);args.output = os.path.abspath(args.output)



    python_script = os.path.join(_root_dir,'src','ATAC','1.preprocess','data.py')

    cmd = [_CDesk_env['python3_7'], python_script,'-i',args.input,'-o',args.output]

    subprocess.run(cmd,check=True,env=_CDesk_env)



    bash_script = os.path.join(_root_dir,'src','ATAC','1.preprocess','ATACseqData.sh')

    cmd = [bash_script,'-i',os.path.join(args.output,'data'),'-o',args.output,'-s',args.species,'-p',str(args.thread),'-g',args.input]

    subprocess.run(cmd,check=True,env=_CDesk_env)

# ATAC: 2. QC
def ATAC_cor(args):
    python_script = os.path.join(_root_dir,'src','ATAC','2.correlation','cor.py')
    args.input = args.input.rstrip('/');args.output = args.output.rstrip('/')
    args.input = os.path.abspath(args.input);args.output = os.path.abspath(args.output)
    cmd = [_CDesk_env['python3_7'], python_script,'-i',args.input,'--species',args.species,'--bin',str(args.bin),
           '--step',str(args.step),'-o',args.output,'-t',str(args.thread),'--width',str(args.width),'--height',str(args.height)]
    subprocess.run(cmd,check=True,env=_CDesk_env)

# ATAC: 3. Dynamic
def ATAC_dynamic(args):
    python_script = os.path.join(_root_dir,'src','ATAC','3.dynamic','dynamic.py')
    args.input = args.input.rstrip('/');args.output = args.output.rstrip('/')
    args.input = os.path.abspath(args.input);args.output = os.path.abspath(args.output)
    cmd = [_CDesk_env['python3_7'], python_script,'--file_meta',args.input,
           '--remove_cluster',args.remove_cluster,'--important_index',str(args.important_index),
           '--species',args.species,'-o',args.output]
    subprocess.run(cmd,check=True,env=_CDesk_env)

# ATAC 4. motif
def ATAC_motif(args):
    if args.motif_mode =='homer':
        python_script = os.path.join(_root_dir,'src','ATAC','4.motif','FindMotif_homer.py')
        args.input = args.input.rstrip('/');args.output = args.output.rstrip('/')
        args.input = os.path.abspath(args.input);args.output = os.path.abspath(args.output)
        cmd = [_CDesk_env['python3_7'], python_script,args.input,args.output,args.species,str(args.thread),args.mode,str(args.start),str(args.end)]
        subprocess.run(cmd,check=True,env=_CDesk_env)

    elif args.motif_mode =='meme':
        python_script = os.path.join(_root_dir,'src','ATAC','4.motif','FindMotif_meme.py')
        args.input = args.input.rstrip('/');args.output = args.output.rstrip('/');args.database=args.database.rstrip('/')
        args.input = os.path.abspath(args.input);args.output = os.path.abspath(args.output);args.database=os.path.abspath(args.database)
        cmd = [_CDesk_env['python3_7'], python_script,args.input,args.database,args.output,str(args.thread),str(args.nmotifs)]
        subprocess.run(cmd,check=True,env=_CDesk_env)

    elif args.motif_mode =='homer_heatmap':
        R_script = os.path.join(_root_dir,'src','ATAC','4.motif','homer_heatmap.R')
        args.input = args.input.rstrip('/');args.output = args.output.rstrip('/')
        args.input = os.path.abspath(args.input);args.output = os.path.abspath(args.output)
        if args.motif != 'no':
            args.motif = args.motif.rstrip('/');args.motif = os.path.abspath(args.motif)
        cmd = ["Rscript",R_script,args.input,args.output,args.motif,str(args.width),str(args.height)]
        subprocess.run(cmd,check=True,env=_CDesk_env)

    elif args.motif_mode =='enrich_heatmap':
        python_script = os.path.join(_root_dir,'src','ATAC','4.motif','motif_enrich_heatmap.py')
        args.motif_bed = args.motif_bed.rstrip('/');args.output = args.output.rstrip('/');args.bed = args.bed.rstrip('/')
        args.motif_bed = os.path.abspath(args.motif_bed);args.output = os.path.abspath(args.output);args.bed = os.path.abspath(args.bed)
        cmd = [_CDesk_env['python3_7'], python_script,'--motif_bed',args.motif_bed,'-o',args.output,'--bed',args.bed,
               '--species',args.species,"-t",str(args.thread)]

        subprocess.run(cmd,check=True,env=_CDesk_env)

# ATAC: 5. Accessbility
def ATAC_accessbility(args):
    python_script = os.path.join(_root_dir,'src','ATAC','5.accessbility','accessbility.py')
    args.input = args.input.rstrip('/');args.output = args.output.rstrip('/');args.bw=args.bw.rstrip('/')
    args.input = os.path.abspath(args.input);args.output = os.path.abspath(args.output);args.bw = os.path.abspath(args.bw)
    if args.region_list != '':
        args.region_list = args.region_list.rstrip('/');args.region_list = os.path.abspath(args.region_list)
    if args.mode == 'gene':
        if args.species == None:
            print("Please provide gene parameter for gene mode")
            sys.exit(1)
        cmd = [_CDesk_env['python3_7'], python_script,'-bw',args.bw,'-o',args.output,'-m',args.mode,
               '-i',args.input,'-s',args.species,'--promoter_size',str(args.promoter_size),
               '--method',args.method,'--region_list',args.region_list,'--center',args.center]
    else:
        cmd = [_CDesk_env['python3_7'], python_script,'-bw',args.bw,'-o',args.output,'-m',args.mode,
               '-i',args.input,'--promoter_size',str(args.promoter_size),
               '--method',args.method,'--region_list',args.region_list,'--center',args.center]
    subprocess.run(cmd,check=True,env=_CDesk_env)

# ATAC: 6. IGV snapshots
def ATAC_igv(args):
    python_script = os.path.join(_root_dir,'src','ATAC','6.igv','snapshot.py')
    args.input = args.input.rstrip('/');args.output = args.output.rstrip('/');args.region=args.region.rstrip('/')
    args.input = os.path.abspath(args.input);args.output = os.path.abspath(args.output);args.region=os.path.abspath(args.region)
    cmd = [_CDesk_env['python3_6'], python_script,'-i',args.input,'-o',args.output,'--genome',args.genome,
               '--display_mode',args.display_mode,'--region',args.region,
               '--type',args.type]
    subprocess.run(cmd,check=True,env=_CDesk_env)

# ATAC: 7. Similarity
def ATAC_similarity(args):
    bash_script = os.path.join(_root_dir,'src','ATAC','7.similarity','ATAC_data_correlation.sh')
    args.bw1 = args.bw1.rstrip('/');args.output = args.output.rstrip('/');args.bw2 = args.bw2.rstrip('/')
    args.bw1 = os.path.abspath(args.bw1);args.output = os.path.abspath(args.output);args.bw2 = os.path.abspath(args.bw2)
    if args.region:
        args.region = args.region.rstrip('/');args.region = os.path.abspath(args.region)
        cmd = [bash_script,args.bw1,args.bw2,
            args.output,str(args.thread),args.region]
    else:
        cmd = [bash_script,args.bw1,args.bw2,
            args.output,str(args.thread)]
    subprocess.run(cmd,check=True,env=_CDesk_env)

def handle(args):
    """ATAC command dispatcher"""
    cmd = args.ATAC_command
    if cmd == 'preprocess':
        ATAC_preprocess(args)
    elif cmd == 'correlation':
        ATAC_cor(args)
    elif cmd == 'dynamic':
        ATAC_dynamic(args)
    elif cmd == 'motif':
        ATAC_motif(args)
    elif cmd == 'accessbility':
        ATAC_accessbility(args)
    elif cmd == 'igv':
        ATAC_igv(args)
    elif cmd == 'similarity':
        ATAC_similarity(args)