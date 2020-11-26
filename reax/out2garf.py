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

description=r"""
Author: Jiayan Xu, jxu15@qub.ac.uk
"""

def read_outcar(outcar='OUTCAR', natoms=100 , nframes=1000, wdat=False):
    # check file existence
    if not os.path.exists(outcar):
        raise ValueError('%s doesnot exist.' %outcar)

    # read OUTCAR
    frames = []
    energies = []
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
        if line.startswith('  FREE ENERGIE'):
            fopen.readline() # segment line ---...---
            fopen.readline() # free energy TOTEN
            fopen.readline() # blank line
            data = fopen.readline()
            energies.append(data.strip().split()[-1])
            count += 1
        if count == nframes:
            flag = False
    fopen.close()
    print('Successfully read %s, get positions, forces, energies ...' %outcar)

    energies = np.array(energies, dtype=float)
    
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
    return frames, energies

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

def write_biogrf(name, lx, ly, lz, alpha, beta, gamma, atoms, coords):
    """Generate a lammps-like bgf file."""
    # first few lines
    content = 'BIOGRF 200\n' # FILTYP (A6,I5)
    content += 'DESCRP ' + name + '\n' # ('DESCRP',1X,A8)
    content += 'REMARK generated by pos2bgf.py\n' # ('REMARK',1X,A)

    # ('FORCEFIELD',1X,A8)
    # ('PERIOD',1X,3I1)

    content += 'AXES   xyz\n' # ('AXES',3X,A)
    content += ('CRYSTX ' + 6*'{:>11.5f}' + '\n')\
            .format(lx, ly, lz, alpha, beta, gamma) # ('CRYSTX',1X,6F11.5)
    
    # constraint
    #content += 'BOND RESTRAINT 1 2 0.50 7500.00 7500.00 0.0000000       0       0\n'
    
    # ('ATOM'|'HETATM',1X,I5,1X,A5,1X,A3,1X,A1,1X,A5,3F10.5,1X,A5,I3,I2,1X,F8.5)  
    atom_format = 'HETATM {:>5d} {:>5s} {:>3s} {:>1s} {:>5s}'+3*'{:>10.5f}'+\
            ' {:>5s}{:>3d}{:>2d} {:>8.5f}\n'
    for i, (atom, coord) in enumerate(zip(atoms, coords)):
        content += atom_format.format(i+1, atom, '', '', '', \
                *coord, atom, 0, 0, 0)

    # END
    content += 'END\n'

    return content

def calc_trans(lattice):
    # lattice
    a, b, c = lattice
    lx, ly, lz = norm(a), norm(b), norm(c) # angstrom

    alpha = arccos(dot(b,c)/ly/lz)/pi*180
    beta = arccos(dot(a,c)/lx/lz)/pi*180
    gamma = arccos(dot(a,b)/lx/ly)/pi*180 # degree

    vol = dot(a, cross(b, c))

    lat_lmp = np.array([[lx, 0, 0,], \
            [ly*np.cos(gamma/180*np.pi), ly*np.sin(gamma/180*np.pi), 0], \
            [0, 0, lz]])

    a_lmp, b_lmp, c_lmp = lat_lmp[0], lat_lmp[1], lat_lmp[2]
    trans_matrix = 1/vol*np.dot(lat_lmp.T, \
            [np.cross(b,c),np.cross(c,a),np.cross(a,b)])
    return lx, ly, lz, alpha, beta, gamma, trans_matrix

def adjust_poses(poses, refposes):
    """poses, poses_ref should be in direct"""
    # atoms, take last step as reference
    for i in range(len(poses)):
        for x in range(3):
            move = round(poses[i][x] - refposes[i][x], 0)
            poses[i][x] -= move
    refposes = poses.copy()

    return poses, refposes

def write_trainsetin(names, energies, refname):
    EVTOKCAM = 23.061
    content = 'ENERGY 100.0\n'
    content += '# weight / structure / reference / DFT in kcal/mol\n'
    sample_format = '1.0   + {:>12s}/1 - {:>12s}/1 {:>12.4f}\n'
    for name, energy in zip(names, energies):
        content += sample_format.format(name, refname, energy*EVTOKCAM)
    content += 'ENDENERGY'

    return content

def out2garf(outcar='OUTCAR', poscar='POSCAR', nframes=100, intv=5,\
        refstructure=['POSCAR', '0.0'], samplename='POSCAR'):
    """
    """
    # read POSCAR
    fname, scaling, lattice, \
            symbols, numbers, atoms, refpos, fixes = read_poscar(poscar)
    natoms = np.sum(numbers)

    lx, ly, lz, alpha, beta, gamma, trans_matrix = calc_trans(lattice)
    
    # read OUTCAR
    frames, energies = read_outcar(outcar, natoms, nframes)

    # adjust positiosn and get cartesian coordinates
    dirposes = []
    cartposes = []
    for frame in frames:
        dirpos = dot(frame[0], inv(lattice.T))
        dirpos, refpos = adjust_poses(dirpos, refpos)
        dirposes.append(dirpos)
        cartposes.append(dot(dirpos, lattice))
    cartposes = np.array(cartposes)

    # write geo
    names = []
    geo_content = ''
    for i, cartpos in enumerate(cartposes):
        if i%intv == 0:
            name = samplename + '_' + str(i+1).zfill(4)
            coords = (dot(trans_matrix, cartpos.T)).T
            geo_content += write_biogrf(name, lx, ly, lz, \
                alpha, beta, gamma, atoms, coords) + '\n'
            names.append(name)

    refposcar, refenergy = refstructure[0], float(refstructure[1])
    if not (isinstance(refposcar, str) and isinstance(refenergy, float)):
        raise ValueError('First must be POSCAR path and second must be energy.')
    fname, scaling, lattice, \
            symbols, numbers, atoms, refpos, fixes = read_poscar(refposcar)
    lx, ly, lz, alpha, beta, gamma, trans_matrix = calc_trans(lattice)

    refname = samplename + '_REF'
    cartpos = dot(refpos, lattice)
    coords = (dot(trans_matrix, cartpos.T)).T
    ref_content = write_biogrf(refname, lx, ly, lz, \
            alpha, beta, gamma, atoms, coords)

    geo_content = ref_content + '\n' + geo_content

    with open('geo', 'w') as writer:
        writer.write(geo_content)
    print('Successfully write geo ...')

    # write trainset.in file
    energies = energies - refenergy
    tsi_content = write_trainsetin(names, energies, refname)

    with open('trainset.in', 'w') as writer:
        writer.write(tsi_content)
    print('Successfully write trainset.in ...')


if __name__ == '__main__':
    # args
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-v', '--version', \
            action='version', version='%(prog)s 0.1')
    parser.add_argument('-o', '--outcar', nargs='?',\
            default='OUTCAR', help='OUTCAR')
    parser.add_argument('-p', '--poscar', nargs='?',\
            default='POSCAR', help='POSCAR')
    parser.add_argument('-nf', '--nframes', required=True,\
            type=int, help='Number of Frames')
    parser.add_argument('-i', '--interval', \
            type=int, default=10, help='Selection Interval')
    parser.add_argument('-r', '--refstr', nargs=2, required=True,\
            help='Reference Structure')
    parser.add_argument('-sn', '--samplename', required=True,\
            help='Sample Name')
    args = parser.parse_args()
    
    #out2garf(outcar='OUTCAR', poscar='POSCAR', nframes=5000, intv=500,\
    #        refstructure=['POSCAR', '-5000.0'], samplename='CO2ad')
    out2garf(args.outcar, args.poscar, args.nframes, args.interval,\
            args.refstr, args.samplename)
