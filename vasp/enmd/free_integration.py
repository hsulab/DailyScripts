#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import time
import argparse

import numpy as np
from scipy import integrate

pi = np.pi

norm = np.linalg.norm
inv = np.linalg.inv

dot = np.dot
cross = np.cross

arccos = np.arccos

description=r"""
Jiayan Xu, jxu15@qub.ac.uk - 
    This script for work integration from free MD.
Bug Fix:
    2020.2 Adjust positions by PBC.
"""


def read_outcar(outcar='OUTCAR', natoms=100 , nframes=1000, wdat=False):
    # read OUTCAR
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
            frames.append([poses, forces])
            count += 1
        if count == nframes:
            flag = False
    fopen.close()
    print('Successfully read OUTCAR, get positions and forces ...')
    
    # out datfile ?
    '''
    for i, d in enumerate(data):
        content = '#\n'
        for j in range(nsteps):
            pos, force = data[i][0][j], data[i][1][j]
            content += ('{:<12.4f}'*6+'\n').format(*pos, *force)
        with open('atom-'+str(i+1)+'.dat', 'w') as writer:
            writer.write(content)
        break
    '''
    return frames

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

    print('Successfully read POSCAR, taking it as the reference...')
    return fname, scaling, lattice, symbols, numbers, poses, fixes

def adjust_poses(poses, refposes):
    """poses, poses_ref should be in direct"""
    # atoms, take last step as reference
    for i in range(len(poses)):
        for x in range(3):
            move = poses[i][x] - refposes[i][x]
            move = round(move, 0)
            poses[i][x] -= move
    refposes = poses.copy()

    return poses, refposes



def atom_integration(poses, forces, nsteps, fix):
    e_atom = []
    for i in range(nsteps):
        cur_fracen = []
        for x in range(3):
            if fix[x] == 'T':
                pos, force = poses[:,x], forces[:,x]
                # work and energy opposite
                v = -integrate.trapz(force[:i+1], pos[:i+1])
            elif fix[x] == 'F':
                v = 0.0
            else:
                raise ValueError('Wrong fix')
            cur_fracen.append(v)
        e_sum = np.sum(cur_fracen)
        cur_fracen.append(e_sum)
        e_atom.append(cur_fracen)
    e_atom = np.array(e_atom)
    return e_atom

def free_integration(poscar='POSCAR', outcar='OUTCAR', nframes=1):
    # read poscar and outcar
    fname, scaling, lattice, symbols, numbers, refposes, fixes = read_poscar(poscar)
    natoms = np.sum(numbers)

    frames = read_outcar(outcar=outcar, natoms=natoms, nframes=nframes)

    for i, frame in enumerate(frames):
        # adjust positions
        poses = frame[0] # cartesian coordinates
        dirposes = dot(poses, inv(lattice))
        dirposes, refposes = adjust_poses(dirposes, refposes)
        poses = dot(dirposes, lattice)
        frames[i][0] = poses

    # rearange data into atom
    data = [[[],[]] for i in range(natoms)] # atoms positions and forces separate
    for frame in frames:
        poses, forces = frame
        for i in range(natoms):
            data[i][0].append(poses[i])
            data[i][1].append(forces[i])

    for i, d in enumerate(data):
        data[i][0] = np.array(data[i][0])
        data[i][1] = np.array(data[i][1])

    # atoms need integration
    en_atoms = []
    for i, fix in enumerate(fixes):
        poses, forces = data[i]
        en_atom = atom_integration(poses, forces, nframes, fix)
        en_atoms.append(en_atom)
        print('atom %d' %(i+1))
    
    # atoms
    en_sum = np.zeros((nframes,4))
    atom_content = '# Free Molecular Dynamics Integration for Atoms\n'
    for i, (en_atom, fix) in enumerate(zip(en_atoms, fixes)):
        atom_content += ('{:<8s}'+'{:>12s}'*3+'\n')\
                .format('# %6d' %(i+1), *fix)
        atom_content += ('{:<8s}'+'{:>12s}'*4+'\n')\
                .format('  # Step', 'Ex', 'Ey', 'Ez', 'E')
        for j, en in enumerate(en_atom):
            atom_content += ('{:>8d}'+'{:>12.4f}'*4+'\n').format(j+1, *en)
        atom_content += '\n'
        en_sum += en_atom

    with open('FI-ATOM.dat', 'w') as writer:
        writer.write(atom_content)

    # summation
    content = '# Free Molecular Dynamics Integration Summation\n'
    content += ('# SUM\n')
    content += ('{:<8s}'+'{:>12s}'*4+'\n')\
        .format('  # Step', 'Ex', 'Ey', 'Ez', 'E')
    for i, en in enumerate(en_sum):
        content += ('{:>8d}'+'{:>12.4f}'*4+'\n').format(i+1, *en)
    content += '\n'

    with open('FI-SUM.dat', 'w') as writer:
        writer.write(content)

    print('Successfully writer integration data in FI.dat.')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.1')
    parser.add_argument('-p', '--pos', nargs='?', default='POSCAR', help='POSCAR File')
    parser.add_argument('-o', '--out', nargs='?', default='OUTCAR', help='OUTCAR File')
    parser.add_argument('-nf', '--nframes', nargs='?', default=2000,\
            required=True, type=int, help='Number of Frames')

    args = parser.parse_args()

    free_integration(args.pos, args.out, args.nframes)
    #free_integration('POSCAR.ref', 'OUTCAR', 2000)
