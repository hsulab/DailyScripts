#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys

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
                nframes = int(kwargs[nframes])

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
            frames.append((poses, forces))
            count += 1
        if count == nframes or (not line):
            flag = False
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

    return frames

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

def write_arc(frames,refposes,atoms,lattice,cell_para,trans,arcname):
    content = '!BIOSYM archive 3\nPBC=ON\n'
    for i, frame in enumerate(frames):
        # adjust positions
        poses = frame[0] # cartesian coordinates
        dirposes = dot(poses, inv(lattice))
        dirposes, refposes = adjust_poses(dirposes, refposes)
        poses = dot(dirposes, lattice)

        # write time
        content += ('%80.4f\n' % 0.0)
        content += '!DATE     %s\n' \
            %time.asctime( time.localtime(time.time()))

        # write crystal
        content += ('PBC'+'{:>10.4f}'*6+'\n').format(*cell_para)

        for atom, coord in zip(atoms, poses):
            coord = dot(trans, coord.T)
            content += '%2s%16.9f%16.9f%16.9f%5s%2d%8s%8s%7.3f\n' %\
                (atom, coord[0], coord[1], coord[2],
                'XXXX', 1, 'xx', atom, 0.0)
        content += 'end\nend\n'
        print('Writing frame %d ...' %(i+1), end='\r')
    print('Writing frame %d ...' %(i+1))

    arcname = os.path.basename(os.getcwd()) + '-ref.arc'
    with open(arcname, 'w') as writer:
        writer.write(content)

    print('Successfully write positions in OUTCAR to %s.' %arcname)

    return

def out2arc(outcar='OUTCAR',poscar='POSCAR',nframes=100,arcname='vasp.arc'):
    """ outcar to arc"""
    # read POSCAR
    fname, scaling, lattice, symbols, numbers, refposes, fixes = \
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

    # symbols
    atoms = []
    for s, n in zip(symbols, numbers):
        atoms.extend([s]*n)

    # read OUTCAR and write
    frames = read_outcar(outcar=outcar,natoms=natoms,nframes=nframes)

    write_arc(frames,refposes,atoms,lattice,cell_para,trans_matrix,arcname)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-v', '--version', action='version', \
            version='%(prog)s 1.2')
    parser.add_argument('-p', '--pos', nargs='?', default='POSCAR', \
            help='POSCAR File')
    parser.add_argument('-o', '--out', nargs='?', default='OUTCAR', \
            help='OUTCAR File')
    parser.add_argument('-a', '--arc', nargs='?', default='vasp.arc', \
            help='ARC File')
    parser.add_argument('-nf', '--nframes', nargs='?', default=-100, \
            type=int, help='Number of Frames')

    args = parser.parse_args()

    out2arc(args.out,args.pos,args.nframes)

