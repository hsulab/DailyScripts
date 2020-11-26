#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse

import numpy as np

from ase import Atoms
from ase.io import read, write

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-e', '--element', 
        nargs = 2, type=str, 
        required = True, 
        help='elements for the dimer'
    )
    parser.add_argument(
        '-d', '--distance', 
        nargs = 3, type=float, 
        default = [0.8,2.0,13], 
        help='dimer distance'
    )

    parser.add_argument(
        '-o', '--output', 
        default = 'dimer.xyz', 
        help='dimer output xyz'
    )

    args = parser.parse_args()

    elements = args.element

    cell = np.array(
        [[10.,0.,0.],[0.,10.,0.],[0.,0.,10.]]
    )

    dis_range = args.distance
    dis_range[2] = int(dis_range[2])
    distances = np.linspace(*dis_range)
    print(distances)

    frames = []
    for d in distances:
        atoms = Atoms(
            ''.join(elements), 
            positions = [[0.,0.,0.],[0.,0.,d]],
            pbc = [1,1,1],
            cell = cell, 
        )
        frames.append(atoms)

    write(args.output, frames)
