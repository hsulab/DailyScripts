#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from pathlib import Path

import numpy as np
from ase.io import read, write

"""
bulk energy
surface unrelaxed and relaxed energies
surface area
"""

parser = argparse.ArgumentParser()
parser.add_argument(
    '-d', '--bulk', 
    default='./bulk', help="bulk calculation directory"
)
parser.add_argument(
    '-s', '--surf', 
    default='./surf', help="surface calculation directory"
)

args = parser.parse_args()

frames = read(Path(args.bulk)/'OUTCAR', ':')
bulk_energy = frames[-1].get_potential_energy()
natoms_bulk = len(frames[-1])
print("natoms in bulk: ", natoms_bulk)

frames = read(Path(args.surf)/'OUTCAR', ':')
unrelaxed_energy = frames[0].get_potential_energy()
relaxed_energy = frames[-1].get_potential_energy()
cell = frames[-1].cell
surface_area = np.linalg.norm(np.cross(cell[0], cell[1]))
natoms_surface = len(frames[-1])
print("natoms in surface: ", natoms_surface)

coef = natoms_surface / natoms_bulk
unrelaxed_surfen = 0.5 * (unrelaxed_energy - coef*bulk_energy) / surface_area
relaxed_surfen = (0.5 * (unrelaxed_energy - coef*bulk_energy) + (relaxed_energy-unrelaxed_energy)) / surface_area

print('Surface ', args.bulk, args.surf)
print("Surface Area [AA^3]: ", surface_area)
eV2J = 16.02
print('surface energy ([eV/AA2], [J/m2]): %.4f %.4f' %(relaxed_surfen, relaxed_surfen*eV2J))
print('unrelaxed surface energy ([eV/AA2], [J/m2]): %.4f %.4f' %(unrelaxed_surfen, unrelaxed_surfen*eV2J))


if __name__ == '__main__':
    pass
