#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import time
import argparse
from pathlib import Path

from tqdm import tqdm
from ase.io import read, write

from joblib import Parallel, delayed

#d = Path('./PtO_200_500')
#d = Path('./bPtO2_200_500')

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
    '-nj', '--njobs', type=int,
    default=4, help='upper limit on number of directories'
)
parser.add_argument(
    '-l', '--limit', type=int,
    default=10000, help='upper limit on number of directories'
)

args = parser.parse_args()

def find_vasp_dirs(wd):
    cur_vasp_dirs = []
    for p in wd.glob(args.pattern):
        cur_vasp_dirs.append(p)
    print('find number of vasp dirs %d in %s' %(len(cur_vasp_dirs), wd))

    return cur_vasp_dirs

d = Path(args.dir)

vasp_dirs = []
for p in d.parent.glob(d.name+'*'):
    if p.is_dir():
        vasp_dirs.extend(find_vasp_dirs(p))
print('total vasp dirs: %d' %(len(vasp_dirs)))

vasp_dirs_sorted = sorted(vasp_dirs, key=lambda k: int(k.name.split('_')[-1])) # sort by name
print(vasp_dirs_sorted)

def extract_atoms(p):
    vasprun = Path(p) / args.vaspfile
    #atoms = read(vasprun, format='vasp-xml')
    atoms = read(vasprun)

    return atoms

st = time.time()

if args.njobs > 1:
    print('using num of jobs: ', args.njobs)
    frames = Parallel(n_jobs=args.njobs)(delayed(extract_atoms)(p) for p in vasp_dirs_sorted)
else:
    frames = []
    for idx, p in tqdm(enumerate(vasp_dirs_sorted)):
        if idx >= args.limit:
            break
        vasprun = Path(p) / args.vaspfile
        #atoms = read(vasprun, format='vasp-xml')
        atoms = read(vasprun)
        frames.append(atoms)

et = time.time()
print(et-st)

write(d.name+'_sorted.xyz', frames)

if __name__ == '__main__':
    pass
