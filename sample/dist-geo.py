#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

norm = np.linalg.norm

def read_gen(gen='geo_in.gen'):
    with open(gen,'r') as reader:
        lines = reader.readlines()

    # first line
    data = lines[0].strip().split()
    natoms, coordtype = int(data[0]), data[1]

    # second line, symbols
    data = lines[1].strip().split()
    symbols = data

    # poses
    data = [line.strip().split() for line in lines[2:2+natoms]]
    poses = np.array(data,dtype=float)[:,2:]

    # lattice
    data = [line.strip().split() for line in lines[2+natoms+1:2+natoms+4]]
    lattice = np.array(data,dtype=float)

    return natoms, coordtype, symbols, poses, lattice


if __name__ == '__main__':
    # read gen file
    natoms, coordtype, symbols, poses, lattice = read_gen('geo_in.gen')

    # calc dist
    vec = poses[0] - poses[1]
    cur_dis = norm(np.dot(lattice.T,vec))

    print(cur_dis)
