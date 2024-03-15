#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from pathlib import Path
import numpy as np

from ase.io import read, write
from ase.geometry import find_mic

parser = argparse.ArgumentParser()
parser.add_argument(
    "-f", "--file", default="POSCAR",
    help = "vasp file name"
)
args = parser.parse_args()

neb_dir = Path("./")

dir_paths = []
for p in neb_dir.glob("[0-9][0-9]"):
    dir_paths.append(p)
dir_paths.sort()

frames = []
for p in dir_paths:
    atoms = read(p / args.file)
    frames.append(atoms)
if args.file == "CONTCAR":
    write("traj_opt.xyz", frames)
elif args.file == "POSCAR":
    write("traj-ini.xyz", frames)

nframes = len(frames)
differences = np.zeros(len(frames))
init_pos = frames[0].get_positions()
for i in range(1,nframes):
    a = frames[i]
    vector = a.get_positions() - init_pos
    vmin, vlen = find_mic(vector, a.get_cell())
    differences[i] = np.linalg.norm(vlen)

neb_dat = neb_dir / "neb.dat"
if neb_dat.exists():
    data = np.loadtxt(neb_dir / "neb.dat")
    data[:, 1] = differences
    np.savetxt(neb_dir / "MyNeb.dat", data)


if __name__ == '__main__':
    pass
