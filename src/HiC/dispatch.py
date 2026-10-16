"""
HiC functions
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

# HiC: 1. Preprocess
def HiC_preprocess(args):
    args.input = args.input.rstrip('/');args.output = args.output.rstrip('/')
    args.input = os.path.abspath(args.input);args.output = os.path.abspath(args.output)
    python_script = os.path.join(_root_dir,'src','HiC','1.Preprocess','data.py')
    cmd = [_CDesk_env['python3_7'], python_script,'-i',args.input,'-o',args.output]
    subprocess.run(cmd,check=True,env=_CDesk_env)

    bash_script = os.path.join(_root_dir,'src','HiC','1.Preprocess','HiCData.sh')
    cmd = [bash_script,'-i',os.path.join(args.output,'data'),'-o',args.output,'-s',args.species,'-t',str(args.thread),'-f',str(args.filter),'-r',args.resolution]
    subprocess.run(cmd,check=True,env=_CDesk_env)

# HiC: 2. Correlation
def HiC_cor(args):

    python_script = os.path.join(_root_dir,'src','HiC','2.Correlation','correlation.py')

    args.input = args.input.rstrip('/');args.output = args.output.rstrip('/')

    args.input = os.path.abspath(args.input);args.output = os.path.abspath(args.output)

    cmd = [_CDesk_env['python3_7'], python_script,'-i',args.input,'-o',args.output,'--chr',args.chr,

            '--inter_chrom',str(args.inter_chrom),'--width',str(args.width),'--height',str(args.height)]

    subprocess.run(cmd,check=True,env=_CDesk_env)

# HiC: 3. Matrix
def HiC_matrix(args):

    args.input = args.input.rstrip('/');args.output = args.output.rstrip('/')

    args.input = os.path.abspath(args.input);args.output = os.path.abspath(args.output)

    python_script = os.path.join(_root_dir,'src','HiC','3.Matrix','matrix.py')

    cmd = [_CDesk_env['python3_7'], python_script,'-i',args.input,'-o',args.output,'--inputFormat',args.inputFormat,

               '--outputFormat',args.outputFormat,'-t',str(args.thread),'-n',args.norm,

               '-r',args.resolution,'--oe',args.oe,'--sparse',args.sparse,'--species',args.species]

    subprocess.run(cmd,check=True,env=_CDesk_env)

# HiC: 4.TAD
def HiC_tad(args):

    python_script = os.path.join(_root_dir,'src','HiC','4.TAD','TAD.py')

    args.input = args.input.rstrip('/');args.output = args.output.rstrip('/')

    args.input = os.path.abspath(args.input);args.output = os.path.abspath(args.output)

    cmd = [_CDesk_env['python3_7'], python_script,'-i',args.input,'-o',args.output,

               "--insulation_window",str(args.insulation_window),"--directionality_window",str(args.directionality_window),"--expand",str(args.expand)]

    subprocess.run(cmd,check=True,env=_CDesk_env)

# HiC: 5. Compartment
def HiC_compartment(args):
    python_script = os.path.join(_root_dir,'src','HiC','5.Compartment','Compartment.py')
    args.input = args.input.rstrip('/');args.output = args.output.rstrip('/');args.GC = args.GC.rstrip('/')
    args.input = os.path.abspath(args.input);args.output = os.path.abspath(args.output);args.GC = os.path.abspath(args.GC) if args.GC !='no' else 'no'
    cmd = [_CDesk_env['python3_7'], python_script,'-i',args.input,'-r',str(args.resolution),'-o',args.output,
               '-t',str(args.thread),'--GC',args.GC,'--species',args.species,
               '--width',str(args.width),'--height',str(args.height)]
    subprocess.run(cmd,check=True,env=_CDesk_env)

# HiC: 6. Loop
def HiC_loop(args):
    python_script = os.path.join(_root_dir,'src','HiC','6.Loop','loop.py')
    args.input = args.input.rstrip('/');args.output = args.output.rstrip('/')
    args.input = os.path.abspath(args.input);args.output = os.path.abspath(args.output)
    cmd = [_CDesk_env['python3_7'],python_script,'-i',args.input,'-o',args.output,'-t',str(args.thread)]
    subprocess.run(cmd,check=True,env=_CDesk_env)

# HiC: 7. Reconstruct
def HiC_reconstruct(args):
    R_script = os.path.join(_root_dir,'src','HiC','7.Reconstruction','FLAMINGO.R')
    args.input = args.input.rstrip('/');args.output = args.output.rstrip('/')
    args.input = os.path.abspath(args.input);args.output = os.path.abspath(args.output)
    cmd = ['Rscript', R_script,args.input,args.chr,args.domain_res,args.frag_res,args.norm,str(args.thread),args.output]
    subprocess.run(cmd,check=True,env=_CDesk_env)

# HiC: 8. PSS
def HiC_pss(args):
    python_script = os.path.join(_root_dir,'src','HiC','8.PSS','pss.py')
    args.input = args.input.rstrip('/');args.output = args.output.rstrip('/')
    args.input = os.path.abspath(args.input);args.output = os.path.abspath(args.output)
    cmd = [_CDesk_env['python3_7'],python_script,'-i',args.input,'-o',args.output,
            '--bin',args.bin,'--curve_mean_min',args.curve_mean_min,'--curve_bin_min',args.curve_bin_min,'--curve_x_min',args.curve_x_min,'--curve_x_max',args.curve_x_max,
            '--expcurve_y_bias',args.expcurve_y_bias,'--expcurve_x_min',args.expcurve_x_min,'--expcurve_x_max',args.expcurve_x_max,'--heatmap_x_min',args.heatmap_x_min,
            '--heatmap_x_max',args.heatmap_x_max,'--width',args.width,'--height',args.height]
    subprocess.run(cmd,check=True,env=_CDesk_env)

# HiC: 9. Compare
def HiC_compare(args):
    args.hic1 = args.hic1.rstrip('/');args.hic2 = args.hic2.rstrip('/');args.output = args.output.rstrip('/')
    args.hic1 = os.path.abspath(args.hic1);args.hic2 = os.path.abspath(args.hic2);args.output = os.path.abspath(args.output)
    python_script = os.path.join(_root_dir,'src','HiC','9.Contact_compare','compare.py')
    cmd = [_CDesk_env['python3_7'],python_script,'--hic1',args.hic1,'--hic2',args.hic2,'-o',args.output,'--chr',args.chr,'--A_min',str(args.A_min)]
    subprocess.run(cmd,check=True,env=_CDesk_env)

def handle(args):
    """HiC command dispatcher"""
    cmd = args.HiC_command
    if cmd == 'preprocess':
        HiC_preprocess(args)
    elif cmd == 'correlation':
        HiC_cor(args)
    elif cmd == 'matrix':
        HiC_matrix(args)
    elif cmd == 'TAD':
        HiC_tad(args)
    elif cmd == 'compartment':
        HiC_compartment(args)
    elif cmd == 'loop':
        HiC_loop(args)
    elif cmd == 'reconstruct':
        HiC_reconstruct(args)
    elif cmd == 'pss':
        HiC_pss(args)
    elif cmd == 'compare':
        HiC_compare(args)