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

def read_poscar(poscar='POSCAR', format='vasp5'):
    """read POSCAR"""
    # check file existence
    if not os.path.exists(poscar):
        raise ValueError('%s doesnot exist.' %poscar)

    # read poscar
    with open(poscar, 'r') as reader:
        lines = reader.readlines()
    lines = [line.strip().split() for line in lines]

    fname = ' '.join(lines[0]) # file description
    scaling = float(lines[1][0])
    lattice = np.array(lines[2:5], dtype=float)
    if lines[5][0].isdigit():
        raise ValueError('Please use VASP5 POSCAR format.')
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

    # symbols
    atoms = []
    for s, n in zip(symbols, numbers):
        atoms.extend([s]*n)

    print('Successfully read %s ...' %poscar)

    return fname, scaling, lattice, symbols, numbers, atoms, poses, fixes

def write_poscar(poscar, lattice, symbols, numbers, poses, fixes, vels):
    # write poscar
    content = 'POSCAR Created by Python Script\n'
    content += '{:>12.9f}\n'.format(1.0)
    for lat in lattice:
        content += '  {:>12.8f} {:>12.8f} {:>12.8f}\n'.format(*lat)
    content += ('  '+'{:>4s}'*len(symbols)+'\n').format(*symbols)
    content += ('  '+'{:>4d}'*len(numbers)+'\n').format(*numbers)
    content += 'Selective Dynamics\nDirect\n'
    for pos, fix in zip(poses,fixes):
        content += ('  '+'{:>16.12f}'*3+'{:>4s}'*3+'\n').format(*pos,*fix)
    content += '\n'

    #for vel in vels:
    #    content += ('  '+'{:>16.12f}'*3+'\n').format(*vel)
                                                                           
    with open(poscar, 'w') as writer:
        writer.write(content)

    print('Successfully write %s ...' %poscar)
    
    return

def adjust_poses(poses, refposes):
    """poses, poses_ref should be in direct"""
    # atoms, take last step as reference
    for i in range(len(poses)):
        for x in range(3):
            move = round(poses[i][x] - refposes[i][x], 0)
            poses[i][x] -= move
    refposes = poses.copy()

    return poses, refposes


def lat2mat(lattice):
    # lattice
    a, b, c = lattice
    lx, ly, lz = norm(a), norm(b), norm(c) # angstrom

    vol = dot(a, cross(b, c))

    alpha = arccos(dot(b,c)/ly/lz)/pi*180
    beta = arccos(dot(a,c)/lx/lz)/pi*180
    gamma = arccos(dot(a,b)/lx/ly)/pi*180 # degree

    lat_lmp = np.array([[lx, 0, 0,], \
            [ly*np.cos(gamma/180*np.pi), ly*np.sin(gamma/180*np.pi), 0], \
            [0, 0, lz]])

    a_lmp, b_lmp, c_lmp = lat_lmp[0], lat_lmp[1], lat_lmp[2]
    trans_matrix = 1/vol*np.dot(lat_lmp.T, \
            [np.cross(b,c),np.cross(c,a),np.cross(a,b)])

    cell = [lx,ly,lz,alpha,beta,gamma]

    return cell, trans_matrix


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0')
    parser.add_argument('-r', '--refcar', nargs='?', default='POSCAR', \
            help='POSCAR File')
    parser.add_argument('-p', '--poscar', nargs='?', default='POSCAR', \
            help='POSCAR or CONTCAR File')

    args = parser.parse_args()
    
    # read BM.dat file
    fname, scaling, lattice, \
        symbols, numbers, atoms, refposes, fixes = read_poscar(args.refcar)
    cell, trans = lat2mat(lattice)

    fname, scaling, lattice, \
        symbols, numbers, atoms, poses, fixes = read_poscar(args.poscar)

    poses, refposes = adjust_poses(poses,refposes)

    vels = np.zeros((np.sum(numbers),3))
    write_poscar('POSCAR.new', lattice, symbols, numbers, poses, fixes, vels)
