#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path

import numpy as np

from ase.io import read, write

working_directory = Path('./')

nebdirs = []
for nebdir in working_directory.glob('[0-9][0-9]'):
    if nebdir.is_dir():
        nebdirs.append(nebdir)
nebdirs.sort()

frames = []
for nebdir in nebdirs:
    contcar = nebdir / 'CONTCAR'
    #if contcar.exists():
    #    print(contcar)
    #    atoms = read(contcar)
    #    frames.append(atoms)
    #else:
    #    raise ValueError('contcar not found.')
    outcar = nebdir / 'OUTCAR'
    if outcar.exists():
        print(outcar)
        atoms = read(outcar, '-1')
        print(atoms.get_potential_energy())
        print(np.max(np.fabs(atoms.get_forces())))
        frames.append(atoms)

assert len(frames) == len(nebdirs)
write('neb.xyz', frames)


if __name__ == '__main__':
    pass
