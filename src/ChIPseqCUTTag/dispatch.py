"""
ChIPseqCUTTag functions
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

# ChIPseq&CUTTag: 1.Preprocess
def ChIPseq_CUTTag_preprocess(args):
    args.input = args.input.rstrip('/');args.output = args.output.rstrip('/')
    args.input = os.path.abspath(args.input);args.output = os.path.abspath(args.output)
    python_script = os.path.join(_root_dir,'src','ChIPseqCUTTag','1.preprocess','data.py')
    cmd = [_CDesk_env['python3_7'], python_script,'-i',args.input,'-o',args.output]
    subprocess.run(cmd,check=True,env=_CDesk_env)

    bash_script = os.path.join(_root_dir,'src','ChIPseqCUTTag','1.preprocess','ChIPseqCUTTagData.sh')
    cmd = [bash_script,'-i',os.path.join(args.output,'data'),'-o',args.output,'-s',args.species,'-p',str(args.thread),'-g',args.input]
    subprocess.run(cmd,check=True,env=_CDesk_env)

# ChIPseq&CUTTag: 3. Pattern
def ChIPseq_CUTTag_pattern(args):
    python_script = os.path.join(_root_dir,'src','ChIPseqCUTTag','3.pattern','pattern.py')
    args.bw = args.bw.rstrip('/');args.output = args.output.rstrip('/');args.bed=args.bed.rstrip('/')
    args.bw = os.path.abspath(args.bw);args.output = os.path.abspath(args.output);args.bed=os.path.abspath(args.bed)
    cmd = [_CDesk_env['python3_6'], python_script,'-bw',args.bw,'-bed',args.bed,'-o',args.output,
           '--mode',args.mode,'--region',args.region,
            '--height',str(args.height),'-a',str(args.upstream),'-b',str(args.downstream),'-m',str(args.length)]
    subprocess.run(cmd,check=True,env=_CDesk_env)

def handle(args):
    """ChIPseqCUTTag command dispatcher"""
    cmd = args.ChIPseq_CUTTag_command
    if cmd == 'preprocess':
        ChIPseq_CUTTag_preprocess(args)
    elif cmd == 'correlation':
        from ATAC.dispatch import ATAC_cor
        ATAC_cor(args)
    elif cmd == 'pattern':
        ChIPseq_CUTTag_pattern(args)
    elif cmd == 'motif':
        from ATAC.dispatch import ATAC_motif
        ATAC_motif(args)
    elif cmd == 'accessbility':
        from ATAC.dispatch import ATAC_accessbility
        ATAC_accessbility(args)
    elif cmd == 'similarity':
        from ATAC.dispatch import ATAC_similarity
        ATAC_similarity(args)
    elif cmd == 'igv':
        from ATAC.dispatch import ATAC_igv
        ATAC_igv(args)