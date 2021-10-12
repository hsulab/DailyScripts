#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import time
import argparse
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

d = Path(args.dir)

vasp_dirs = []
for p in d.parent.glob(d.name+'*'):
    if p.is_dir():
        vasp_dirs.extend(find_vasp_dirs(p))
print('total vasp dirs: %d' %(len(vasp_dirs)))

print("sorted by last integer number...")
vasp_dirs_sorted = sorted(vasp_dirs, key=lambda k: int(k.name.split('_')[-1])) # sort by name
#print(vasp_dirs_sorted) 
def extract_atoms(p):
    vasprun = Path(p) / args.vaspfile
    #atoms = read(vasprun, format='vasp-xml')
    atoms = read(vasprun, args.indices)
    #print(atoms)

    return atoms

st = time.time()

if args.njobs > 1:
    print('using num of jobs: ', args.njobs)
    atoms = Parallel(n_jobs=args.njobs)(delayed(extract_atoms)(p) for p in vasp_dirs_sorted)
    frames = []
    if isinstance(atoms, list):
        frames.extend(atoms)
    else:
        frames.append(atoms)
else:
    frames = []
    for idx, p in tqdm(enumerate(vasp_dirs_sorted)):
        if idx >= args.limit:
            break
        vasprun = Path(p) / args.vaspfile
        #atoms = read(vasprun, format='vasp-xml')
        atoms = read(vasprun, args.indices)
        if isinstance(atoms, list):
            frames.extend(atoms)
        else:
            frames.append(atoms)

et = time.time()
print('cost time: ', et-st)


if len(frames) > 0:
    if args.check:
        nconverged = 0
        for atoms in frames:
            indices = [a.index for a in atoms if a.position[0]<1.5]
            cons = FixAtoms(indices = indices)
            atoms.set_constraint(cons)
            forces = atoms.get_forces(apply_constraint=True)
            maxforce = np.max(np.fabs(forces))
            if maxforce < 0.05:
                nconverged += 1
        print("number of converged: ", nconverged)
    print("Number of frames: ", len(frames))
    write(d.name+'_sorted.xyz', frames)
else:
    print("No frames...")

if __name__ == '__main__':
    pass
