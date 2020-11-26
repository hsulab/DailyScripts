#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse

import numpy as np

import ase
import ase.io.vasp as ase_vasp
from ase import Atom, Atoms

MAXFRAME = 10000

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


def read_outcar(outcar='OUTCAR',natoms=100,verbose=True,wdat=False,**kwargs):
    """
    in each frame,
    coordinates, forces, stress will be read
    """
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
    while flag:
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
            #print(stress)
        if line.strip().startswith('VOLUME and BASIS'):
            fopen.readline() # segement line
            fopen.readline() # energy-cutoff
            fopen.readline() # volume of cell
            fopen.readline() # title
            cur_lat = []
            for n in range(3):
                data = fopen.readline().strip().split()
                cur_lat.append(data)
            cur_lat = np.array(cur_lat,dtype=float)[:,0:3]
        if line.strip().startswith('POSITION'):
            fopen.readline() # segment line ---...---
            poses, forces = [], []
            for n in range(natoms):
                data = fopen.readline().strip().split()
                poses.append(data[:3]) # x y z
                forces.append(data[3:]) # fx fy fz
            poses = np.array(poses, dtype=float)
            forces = np.array(forces, dtype=float)
            frames.append((cur_lat,poses, forces, stress))
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
    
    if verbose:
        print('Successfully read OUTCAR, get positions and forces ...')

    return frames, energy, free_energy

def write_xyz(symbols,positions,forces,**kwargs):
    """positions in cartesian (AA) and forces in eV/AA"""
    # check number of atoms
    natoms = len(symbols)

    # implented properties
    comment_properties = ['Lattice','pbc','energy','free_energy',\
            'stress','virial']
    atomic_properties = ['species','pos','forces']

    # check args
    comment_content = ''
    if len(kwargs) > 0:
        for key, value in kwargs.items():
            if key in comment_properties:
                if key in ['energy','free_energy']:
                    # float
                    value = float(value)
                    comment_content += ("{:<s}="+"{:<.4f}"+" ") \
                            .format(key,value)
                elif key in ['Lattice','stress','virial']:
                    # list of float properties
                    value = list(np.array(value,dtype=float).ravel())
                    if len(value) != 9:
                        raise ValueError('Lattice/stress/virial must have 9 components.')
                    comment_content += ("{:<s}="+"\""+"{:<.4f} "*len(value)+"\""+" ") \
                            .format(key,*value)
                elif key in ['pbc']:
                    # list of strings
                    comment_content += ("{:<s}="+"\""+"{:<s} "*len(value)+"\""+" ") \
                            .format(key,*value)
    comment_content += 'Properties=species:S:1:pos:R:3:forces:R:3\n'

    # write content
    content = "{:<d}\n".format(natoms)
    content += comment_content

    for i in range(natoms):
        content += ('{:<4s} '+'{:>12.6f} '*6+'\n')\
            .format(symbols[i],*list(positions[i]),*list(forces[i]))

    return content


def frame2xyz(symbols,outcar,molecule=False):
    frames, energy, free_energy = \
        read_outcar(outcar=outcar,natoms=len(symbols),verbose=False)
    lattice = frames[0][0]
    positions = frames[0][1]
    forces = frames[0][2]
    stress = frames[0][3]

    volume = np.dot(np.cross(lattice[0],lattice[1]),lattice[2])
    virial = -stress*volume

    # check if need stress and virial
    if molecule:
        content = write_xyz(symbols,positions,forces,\
                pbc=['T','T','T'],Lattice=lattice,\
                energy=energy[0],free_energy=free_energy[0])
    else:
        content = write_xyz(symbols,positions,forces,\
                pbc=['T','T','T'],Lattice=lattice,\
                energy=energy[0],free_energy=free_energy[0],\
                stress=stress,virial=virial)

    return content


if __name__ == '__main__':
    # new parser
    parser = argparse.ArgumentParser()
    #parser.add_argument('-s', '--symbols', required=True, \
    #        nargs='*', help='chemical symbols in system')
    parser.add_argument('-m', '--molecule', \
            type=bool, help='no stress/virial in molecule')
    args = parser.parse_args()

    #symbols = ['Pt']*4
    #symbols = args.symbols
    fname, scaling, lattice, formula, numbers, refposes, fixes = \
            read_poscar('./POSCAR')

    symbols = []
    for s, n in zip(formula, numbers):
        symbols.extend([s]*n)

    #content = frame2xyz(symbols)
    outcars = []
    for fname in os.listdir('./POSCARs'):
        if fname.startswith('OUTCAR'):
            outcars.append(fname)
    outcars.sort(key=lambda fname:int(fname.split('_')[-1]))

    content = ''
    for outcar in outcars[1::3]:
        print('Read %s ...' %outcar)
        content += frame2xyz(symbols,os.path.join('./POSCARs',outcar),args.molecule)

    with open('data.xyz','w') as writer:
        writer.write(content)

    print('Writer data to data.xyz')

