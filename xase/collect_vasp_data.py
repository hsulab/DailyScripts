#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from pathlib import Path

from tqdm import tqdm
from ase.io import read, write

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
    '-l', '--limit', type=int,
    default=10000, help='upper limit on number of directories'
)

args = parser.parse_args()

d = Path(args.dir)

vasp_dirs = []
for p in d.glob(args.pattern):
    vasp_dirs.append(str(p))
print('find number of vasp dirs %d' %(len(vasp_dirs)))

vasp_dirs_sorted = sorted(vasp_dirs, key=lambda k: int(k.split('_')[-1]))

frames = []
for idx, p in tqdm(enumerate(vasp_dirs_sorted)):
    if idx >= args.limit:
        break
    vasprun = Path(p) / 'vasprun.xml'
    atoms = read(vasprun, format='vasp-xml')
    frames.append(atoms)
write(d.name+'_sorted.xyz', frames)

if __name__ == '__main__':
    pass
