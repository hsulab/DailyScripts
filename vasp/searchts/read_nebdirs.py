#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from pathlib import Path

import numpy as np

from ase import Atoms
from ase.io import read, write

parser = argparse.ArgumentParser()
parser.add_argument(
    "FILE"
)
parser.add_argument(
    "-i", "--indices", default="-1"
)
args = parser.parse_args()

working_directory = Path('./')

vasp_file = args.FILE

nebdirs = []
for nebdir in working_directory.glob('[0-9][0-9]'):
    if nebdir.is_dir():
        nebdirs.append(nebdir)
nebdirs.sort()
print(nebdirs)

frames = []
for nebdir in nebdirs:
    outcar = nebdir / vasp_file
    if outcar.name in ["OUTCAR", "vasprun.xml"]:
        if outcar.exists():
            print(outcar)
            atoms = read(outcar, args.indices)
            if isinstance(atoms, Atoms):
                print("energy: ", atoms.get_potential_energy())
                print("max force: ", np.max(np.fabs(atoms.get_forces(apply_constraint=True))))
                atoms = [atoms]
            else:
                print(f"nframes: {len(atoms)}")
        else:
            raise ValueError("outcar not found.")
    else:
        atoms = [read(outcar)]
    frames.extend(atoms)

#write("neb-"+vasp_file+".xyz", frames)
if args.indices == "-1":
    assert len(frames) == len(nebdirs)
    write(Path.cwd().name+".xyz", frames)
else:
    write(Path.cwd().name+"_nebtraj.xyz", frames)

if __name__ == '__main__':
    pass
