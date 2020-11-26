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

import ase.units

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

def wrap_frame_results():
    return frame_results

def read_outcar(outcar='OUTCAR',natoms=100,verbose=True,wdat=False,**kwargs):
    """
    in each frame,
    coordinates, forces, stress will be read
    """
    # 
    required_properties = ['positions', 'forces', 'energy', 'free_energy']
    additional_properties = ['stress']

    # how many steps to read
    nframes = MAXFRAME
    if kwargs:
        if 'nframes' in kwargs.keys():
            if kwargs['nframes'] > 0:
                nframes = int(kwargs['nframes'])

    # read OUTCAR
    frames = []
    cur_results = {}
    fopen = open(outcar, 'r')
    count, flag = 0, True
    while True:
        line = fopen.readline()
        if line.strip().startswith('in kB'): # kB means kilobar equals 0.1 GPa
            # XX YY ZZ XY YZ ZX
            stress = -np.array([float(d) for d in line.strip().split()[2:]])
            # vigot notation, XX YY ZZ YZ XZ XY
            stress_vigot = stress[[0,1,2,4,5,3]]*1e-1*ase.units.GPa # in eV/AA^3
            stress = []
            for n in [0,5,4,5,1,3,4,3,2]:
                stress.append(stress_vigot[n])
            stress = np.array(stress).reshape(3,3)
            cur_results['stress'] = stress
            #print(stress)
        if line.startswith(' POSITION'):
            fopen.readline() # segment line ---...---
            poses, forces = [], []
            for n in range(natoms):
                data = fopen.readline().strip().split()
                poses.append(data[:3]) # x y z
                forces.append(data[3:]) # fx fy fz
            poses = np.array(poses, dtype=float)
            forces = np.array(forces, dtype=float)
            cur_results['positions'] = poses
            cur_results['forces'] = forces
        if line.strip().startswith('FREE ENERGIE OF THE ION-ELECTRON SYSTEM'):
            fopen.readline() # segment line ---...---
            # free energy F
            data = fopen.readline().strip().split()
            free_energy = float(data[-2])
            cur_results['free_energy'] = free_energy
            fopen.readline() # blank line
            # energy E0
            data = fopen.readline().strip().split()
            energy = float(data[-1])
            cur_results['energy'] = energy

        # check if read one frame successfully 
        for prop in required_properties:
            if prop not in cur_results:
                break
        else:
            count += 1
            frames.append(cur_results)
            cur_results = {}

        if line:
            if count == nframes:
                flag = False
                break
        else:
            flag = False
            break
            # raise ValueError('position and energy not ')

    fopen.close()
    
    if verbose:
        print('Successfully read OUTCAR, get positions and forces ...')
    
    # out datfile ?
    #if wdat:
    #    for i, d in enumerate(data):
    #        content = '#\n'
    #        for j in range(nsteps):
    #            pos, force = data[i][0][j], data[i][1][j]
    #            content += ('{:<12.4f}'*6+'\n').format(*pos, *force)
    #        with open('atom-'+str(i+1)+'.dat', 'w') as writer:
    #            writer.write(content)
    #        break

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

    # symbols
    symbols = []
    for s, n in zip(formula, numbers):
        symbols.extend([s]*n)

    # read OUTCAR and write
    frames = read_outcar(outcar=outcar,natoms=natoms,nframes=nframes)

    stride = 1
    content = ''
    for i, results in enumerate(frames):
        #positions = adjust_coordinates(lattice,)
        if i%stride == 0:
            # TODO: coordinate system, be careful with forces
            #positions = 
            #print(energy[i])
            #print(results['energy'])
            results.update(symbols=symbols)
            results.update(Lattice=lattice)
            results.update(pbc=['T','T','T'])
            results.update(step=i+1)
            if 'stress' in results.keys():
                results.update(virial=-results['stress']*vol)
            content += write_xyz(**results)

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

