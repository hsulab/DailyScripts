#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import time
import argparse
import subprocess
from pathlib import Path

from tqdm import tqdm
from ase.io import read, write
from ase.constraints import FixAtoms

from joblib import Parallel, delayed

import numpy as np


parser = argparse.ArgumentParser()
parser.add_argument(
    '-d', '--dir', 
    default='./', help='vasp calculation directory'
)
parser.add_argument(
    '-p', '--pattern', 
    default='vasp_0_*', help='vasp directory name pattern'
)
parser.add_argument(
    '-f', '--vaspfile', 
    default='vasprun.xml', help='vasp directory name pattern'
)
parser.add_argument(
    '-i', '--indices', default="-1", 
    help="frame indices to read"
)
parser.add_argument(
    '-nj', '--njobs', type=int,
    default=1, help='upper limit on number of directories'
)
parser.add_argument(
    '-l', '--limit', type=int,
    default=10000, help='upper limit on number of directories'
)
parser.add_argument(
    "--check", action="store_true",
    help="check number of converged configurations"
)

args = parser.parse_args()

def find_vasp_dirs(wd):
    cur_vasp_dirs = []
    for p in wd.glob(args.pattern):
        cur_vasp_dirs.append(p)
    print('find number of vasp dirs %d in %s' %(len(cur_vasp_dirs), wd))

    return cur_vasp_dirs

if __name__ == '__main__':
    # wd
    d = Path(args.dir)

    # find structures
    structures_in = list(d.glob("O*.xyz"))[0]
    print(structures_in)
    frames_in = read(structures_in, ":")
    nframes_in = len(frames_in)
    print("Number of input structures: ", nframes_in)
    
    # find vasp files
    vasp_dirs = []
    for p in d.parent.glob(d.name+'*'):
        if p.is_dir():
            vasp_dirs.extend(find_vasp_dirs(p))
    nvaspdirs = len(vasp_dirs)
    print('total vasp dirs: %d' %(nvaspdirs))
    
    if nframes_in != nvaspdirs:
        print("need to resubmit job...")
    else:
        exit()

    # read job script
    job_script = d / "vasp.slurm"
    with open(job_script, "r") as fopen:
        content = fopen.readlines()
    new_content = content
    for i, line in enumerate(content):
        if line.startswith("./compute_by_ase.py"):
            line_list = line.split()
            line_list[-1] = "\"{}:\"\n".format(nvaspdirs-1)
            new_content[i] = " ".join(line_list)
            print(new_content[i])
    new_content = "".join(new_content)

    with open(job_script, "w") as fopen:
        fopen.write(new_content)

    command = 'sbatch vasp.slurm'
    proc = subprocess.Popen(
        command, shell=True, cwd=d,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        encoding = 'utf-8'
    )
    errorcode = proc.wait(timeout=10) # 10 seconds
    if errorcode:
        raise ValueError('Error in resubmitting job...')

    print(''.join(proc.stdout.readlines()))

