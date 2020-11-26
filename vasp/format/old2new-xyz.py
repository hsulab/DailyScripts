#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os

import numpy as np

MAXLINE = 10000
MAXFRAME = 10000

def read_xyz(xyz,natoms):
    # 
    fopen = open(xyz,'r')

    frames = []
    for i in range(MAXLINE):
        line = fopen.readline()
        if line.strip():
            assert int(line.strip().split()[0]) == natoms
            line = fopen.readline() # comment
            forces, poses = [], []
            for j in range(natoms):
                line = fopen.readline() #
                data = line.strip().split()
                poses.append(data[1:4])
                forces.append(data[4:7])
            poses = np.array(poses,dtype=float)
            forces = np.array(forces,dtype=float)
            frames.append([poses,forces])

        if not line:
            break
    else:
        raise ValueError('Too many lines in %s' %xyz)

    fopen.close()

    return frames

def read_outcar(outcar='OUTCAR',natoms=100,verbose=True,wdat=False,**kwargs):
    # how many steps to read
    nframes = MAXFRAME
    if kwargs:
        if 'nframes' in kwargs.keys():
            if kwargs['nframes'] > 0:
                nframes = int(kwargs['nframes'])

    # read OUTCAR
    energy = []
    free_energy = []

    frames = []
    fopen = open(outcar, 'r')
    count, flag = 0, True
    while flag:
        line = fopen.readline()
        if line.startswith(' POSITION'):
            fopen.readline() # segment line ---...---
            poses, forces = [], []
            for n in range(natoms):
                data = fopen.readline().strip().split()
                poses.append(data[:3]) # x y z
                forces.append(data[3:]) # fx fy fz
            poses = np.array(poses, dtype=float)
            forces = np.array(forces, dtype=float)
            frames.append((poses, forces))
            count += 1
        if line.strip().startswith('FREE ENERGIE OF THE ION-ELECTRON SYSTEM'):
            fopen.readline() # segment line ---...---
            # free energy F
            data = fopen.readline().strip().split()
            free_energy.append(float(data[-2]))
            fopen.readline() # blank line
            # energy E0
            data = fopen.readline().strip().split()
            energy.append(float(data[-1]))

        if count == nframes or (not line):
            flag = False
    fopen.close()

    return frames, energy, free_energy

def read_poscar(poscar='POSCAR',format='vasp5',verbose=True):

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

    # TODO: read velocity

    if verbose:
        print('Successfully read POSCAR, taking it as the reference...')

    return fname, scaling, lattice, symbols, numbers, poses, fixes


def write_xyz(symbols,lattice,positions,forces,energy):
    """
    positions in cartesian (AA) and forces in eV/AA
    energy is always force-consistent energy for training purpose
    """
    natoms = len(symbols)

    content = "{:<d}\n".format(natoms)
    content += ("Lattice=\""+"{:<.2f} "*9+"\""+\
            " Properties=species:S:1:pos:R:3:forces:R:3 pbc=\"T T T\""+\
            " energy={:<12.6f}"+"\n") \
            .format(*list(lattice),energy)

    for i in range(natoms):
        content += ('{:<4s} '+'{:>12.6f} '*6+'\n')\
            .format(symbols[i],*list(positions[i]),*list(forces[i]))

    return content


def find_outcars():
    fname, scaling, lattice, formula, numbers, refposes, fixes = \
            read_poscar('POSCAR')

    symbols = []
    for s, n in zip(formula, numbers):
        symbols.extend([s]*n)

    print(symbols)

    #
    outcars = []

    cur_dir = './POSCARs/'
    for cur_fname in os.listdir(cur_dir):
        if cur_fname.startswith('OUTCAR_'):
            outcars.append(cur_fname)

    outcars.sort(key = lambda fname: int(fname.split('_')[-1]))
    print(outcars)

    outcars = outcars[1::3]
    print(outcars)

    content = ''
    for outcar in outcars:
        print(os.path.abspath(outcar))
        outcar = os.path.join(cur_dir,outcar)
        frames, energy, free_energy = read_outcar(outcar=outcar,natoms=np.sum(numbers))
        positions, forces = frames[0][0], frames[0][1]
        en = free_energy[0]
        content += write_xyz(symbols,lattice.ravel(),positions,forces,en)

    with open('data.xyz', 'w') as writer:
        writer.write(content)

    return

if __name__ == '__main__':
    #frames = read_xyz('bonds.xyz',2)
    #print(frames[0])
    find_outcars()
