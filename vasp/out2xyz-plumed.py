#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys

import time
import argparse

import numpy as np
from scipy import integrate

sys.path.append(os.path.join(os.getenv("HOME"), \
        'repository/DailyScripts/common'))
from coreXYZ import write_xyz

pi = np.pi

norm = np.linalg.norm
inv = np.linalg.inv

dot = np.dot
cross = np.cross

arccos = np.arccos

description="""
Jiayan Xu, jxu15@qub.ac.uk -
This script writes arc format file using positions in OUTCAR. During 
convertion, positions will be adjusted according to the reference POSCAR.
"""

MAXFRAME = 100000 # maximum steps in outcar, this is a very large number for safe

def read_outcar(outcar='OUTCAR',natoms=100,verbose=True,wdat=False,**kwargs):
    # how many steps to read
    nframes = MAXFRAME
    if kwargs:
        if 'nframes' in kwargs.keys():
            if kwargs['nframes'] > 0:
                nframes = int(kwargs['nframes'])

    # read OUTCAR
    energy, free_energy = [], []

    frames = []
    fopen = open(outcar, 'r')
    count, flag = 0, True
    if_en, if_pos = False, False
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
            if_pos = True
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
            if_en = True

        if line:
            if if_en and if_pos:
                if_en, if_pos = False, False
                if count == nframes:
                    flag = False
        else:
            flag = False
            # raise ValueError('position and energy not ')

    fopen.close()
    
    if verbose:
        print('Successfully read OUTCAR, get positions and forces ...')
    
    # out datfile ?
    if wdat:
        for i, d in enumerate(data):
            content = '#\n'
            for j in range(nsteps):
                pos, force = data[i][0][j], data[i][1][j]
                content += ('{:<12.4f}'*6+'\n').format(*pos, *force)
            with open('atom-'+str(i+1)+'.dat', 'w') as writer:
                writer.write(content)
            break

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


def adjust_poses(poses,refposes):
    """poses, poses_ref should be in direct"""
    # atoms, take last step as reference
    for i in range(len(poses)):
        for x in range(3):
            move = round(poses[i][x] - refposes[i][x], 0)
            poses[i][x] -= move
    refposes = poses.copy()

    return poses, refposes

def adjust_coordinates():

    return


def out2xyz(outcar='OUTCAR',poscar='POSCAR',descrp='vasp',\
        nframes=100,xyz_fname=None):
    """ outcar to arc"""
    # read POSCAR
    fname, scaling, lattice, formula, numbers, refposes, fixes = \
            read_poscar(args.pos)
    natoms = np.sum(numbers)

    # lattice
    a, b, c = lattice
    lx, ly, lz = norm(a), norm(b), norm(c) # angstrom

    vol = dot(a, cross(b, c))

    alpha = arccos(dot(b,c)/ly/lz)/pi*180
    beta = arccos(dot(a,c)/lx/lz)/pi*180
    gamma = arccos(dot(a,b)/lx/ly)/pi*180 # degree

    cell_para = (lx,ly,lz,alpha,beta,gamma)

    # for hex, prefer a along x and b in x-y plane
    lat_ax = np.array([[lx, 0, 0,], \
            [ly*np.cos(gamma/180*np.pi), ly*np.sin(gamma/180*np.pi), 0], \
            [0, 0, lz]])

    # a_ax, b_ax, c_ax = lat_ax[0], lat_ax[1], lat_ax[2]
    trans_matrix = 1/vol*np.dot(lat_ax.T, \
            [np.cross(b,c),np.cross(c,a),np.cross(a,b)])

    # make sure coordinates are a for x and b in x-y plane

    # read plumed-free forces
    with open('./ORIGINAL_FORCES', 'r') as reader:
        lines = reader.readlines()

    plumedfree_forces = []
    for i, line in enumerate(lines):
        if i%(natoms+2) == 2:
            force_lines = lines[i:i+natoms]
            forces = [f_line.strip().split()[1:4] for f_line in force_lines]
            forces = np.array(forces, dtype=float)
            plumedfree_forces.append(forces)

    # symbols
    symbols = []
    for s, n in zip(formula, numbers):
        symbols.extend([s]*n)

    # read OUTCAR and write
    frames, energy, free_energy = \
            read_outcar(outcar=outcar,natoms=natoms,nframes=nframes)

    stride = 1
    content = ''
    for i, ((positions,xforces), forces) \
            in enumerate(zip(frames,plumedfree_forces)):
        #positions = adjust_coordinates(lattice,)
        if i%stride == 0:
            # TODO: coordinate system, be careful with forces
            #positions = 
            #print(energy[i])
            content += write_xyz(symbols,positions,\
                    forces=forces,\
                    description=descrp,\
                    step=i+1,\
                    pbc=['T','T','T'],Lattice=lattice,\
                    energy=energy[i],free_energy=free_energy[i])

    if not xyz_fname:
        xyz_fname = os.path.basename(os.getcwd()) + '-ref.xyz'

    with open(xyz_fname, 'w') as writer:
        writer.write(content)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-v', '--version', action='version', \
            version='%(prog)s 1.2')
    parser.add_argument('-p', '--pos', nargs='?', default='POSCAR', \
            help='POSCAR File')
    parser.add_argument('-o', '--out', nargs='?', default='OUTCAR', \
            help='OUTCAR File')
    parser.add_argument('-f', '--fname', nargs='?', default=None, \
            help='XYZ Filename')
    parser.add_argument('-d', '--descrp', default='vasp', \
            help='description')
    parser.add_argument('-nf', '--nframes', nargs='?', default=-100, \
            type=int, help='Number of Frames')

    args = parser.parse_args()

    out2xyz(args.out,args.pos,args.descrp,args.nframes,args.fname)

