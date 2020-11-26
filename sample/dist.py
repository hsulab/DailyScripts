#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

norm = np.linalg.norm

def read_poscar(poscar='POSCAR', format='vasp5'):
    """read POSCAR"""
    with open(poscar, 'r') as reader:
        lines = reader.readlines()
    lines = [line.strip().split() for line in lines]

    fname = ' '.join(lines[0]) # file description
    scaling = float(lines[1][0])
    lattice = np.array(lines[2:5], dtype=float)
    symbols = lines[5]
    numbers = [int(i) for i in lines[6]]
    natoms = np.sum(numbers)

    dyntype = ' '.join(lines[7]) # dynamic type
    coorsys = lines[8] # coordinate system

    poses, fixes = [], []
    for coord in lines[9:9+natoms]:
        poses.append(coord[:3])
        fixes.append(coord[3:])
    poses = np.array(poses, dtype=float)

    # print('Successfully read POSCAR, taking it as the reference...')

    return fname, scaling, lattice, symbols, numbers, poses, fixes


if __name__ == '__main__':
    fname, scaling, lattice, symbols, numbers, poses, fixes = read_poscar('POSCAR')
    vec = poses[0] - poses[1]
    dis = norm(np.dot(lattice.T,vec))
    print(dis)
