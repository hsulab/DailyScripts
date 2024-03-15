#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from pathlib import Path

import numpy as np

from ase.io import read, write

working_directory = Path('./')

nebdirs = []
for nebdir in working_directory.glob('[0-9][0-9]'):
    if nebdir.is_dir():
        nebdirs.append(nebdir)
nebdirs.sort()

traj_frames = []
for nebdir in nebdirs[1:-1]: # skip IS and FS
    outcar = nebdir / "OUTCAR"
    if outcar.exists():
        print(outcar)
        frames = read(outcar, ":")
        for i, atoms in enumerate(frames):
            atoms.info["description"] = nebdir.name + "-" + str(i)
        print("nframes: ", len(frames))
        traj_frames.extend(frames)
    else:
        raise ValueError("outcar not found.")

write("neb-traj.xyz", traj_frames)


if __name__ == '__main__':
    pass
