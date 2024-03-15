#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from pathlib import Path

from ase.io import read, write
from ase.io.vasp import write_vasp 
from ase.constraints import FixAtoms

parser = argparse.ArgumentParser()
parser.add_argument(
    "-t", "--traj", 
    help = "neb trajectory"
)
parser.add_argument(
    "-c", "--constraint", type=int, nargs=2,
    help = "constraint index start and end"
)
args = parser.parse_args()

traj = Path(args.traj)

#images = read('../emt-hcp2fcc_top/neb.xyz', '-9:')
#images = read("/mnt/scratch2/users/40247882/oxides/eann-main/train-all/m23r/validations/neb-brg/neb_opt.xyz", ":")
#images = read("/mnt/scratch2/users/40247882/oxides/eann-main/train-all/m23r/validations/neb-top/neb_opt.xyz", ":")
#images = read("/users/40247882/scratch2/oxides/surfaces/Pt111/surf-22/neb/O2diss/neb-fcc-eann_opt.xyz", ":")
#images = read("/users/40247882/scratch2/oxides/surfaces/Pt111/surf-22/neb/O2diss/neb-hcp-eann_opt.xyz", ":")
#images = read('../emt-hcp2fcc_top/neb.xyz', '-9:')
images = read(traj, ":")
print(images)

constraint = FixAtoms(indices=range(args.constraint[0], args.constraint[1]))

for idx, atoms in enumerate(images):
    # create poscar
    poscar_dir = './%s' %(str(idx).zfill(2))
    poscar_dir = Path(poscar_dir)
    poscar_dir.mkdir(parents=True)
    poscar = poscar_dir / 'POSCAR'

    # write atoms
    atoms.set_constraint(constraint)
    write_vasp(poscar, atoms, direct=True, vasp5=True)

if __name__ == '__main__':
    pass
