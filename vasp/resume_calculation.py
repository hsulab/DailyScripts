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
    '-dp', '--dir_pattern', 
    default="O*", help='calculation directory pattern'
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

def find_job_in_queue(d):
    """"""
    command = "squeue --user 40247882 --format=\"%.12i %.12P %.24j %.4t %.12M %.12L %.5D %.4C\""
    proc = subprocess.Popen(
        command, shell=True, cwd=d,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        encoding = 'utf-8'
    )
    errorcode = proc.wait(timeout=10) # 10 seconds
    if errorcode:
        raise ValueError('Error in checking job...')

    lines = proc.stdout.readlines()
    job_ids = [
        line.split()[0] for line in lines
    ][1:]
    #print(job_ids)

    wds = []
    for job_id in job_ids:
        command = "scontrol show job {}".format(job_id)
        proc = subprocess.Popen(
            command, shell=True, cwd=d,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            encoding = 'utf-8'
        )
        errorcode = proc.wait(timeout=10) # 10 seconds
        if errorcode:
            raise ValueError('Error in checking job...')
        lines = proc.stdout.readlines()
        for line in lines:
            line = line.strip()
            if line.startswith("WorkDir="):
                wds.append(line[8:])

    flag = False
    if str(d) in wds:
        print("job already in queue...")
        flag = True

    return flag

def find_vasp_dirs(wd):
    cur_vasp_dirs = []
    for p in wd.glob(args.pattern):
        cur_vasp_dirs.append(p)
    print('find number of vasp dirs %d in %s' %(len(cur_vasp_dirs), wd))

    return cur_vasp_dirs

def check_single_dir(d):
    # wd
    d = Path(d).absolute()

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
        return

    # check job if already in queue
    if find_job_in_queue(d):
        return

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

if __name__ == '__main__':
    cwd = Path(args.dir)
    calc_dirs = list(cwd.glob(args.dir_pattern))
    calc_dirs.sort()
    print("number of dirs: ", len(calc_dirs))
    for d in calc_dirs:
        print("\n\n===== {} =====".format(d))
        check_single_dir(d)
