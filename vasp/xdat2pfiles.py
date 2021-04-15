#!/usr/bin/env python3
from pathlib import Path
import argparse
from tqdm import tqdm
from ase.io import read, write
from ase.io.vasp import read_vasp_xdatcar, write_vasp

parser = argparse.ArgumentParser()
parser.add_argument(
    '-i', '--index', 
    default=':', help='frame index starting from 0 e.g. 0:100'
)

args = parser.parse_args()

frames = read('XDATCAR', index=args.index, format='vasp-xdatcar')

pfiles_path = Path('./pfiles')
pfiles_path.mkdir()
for idx, atoms in tqdm(enumerate(frames)):
    p_path = pfiles_path / ('p'+str(idx+1).zfill(4))
    write_vasp(p_path, atoms, direct=True, vasp5=True)
