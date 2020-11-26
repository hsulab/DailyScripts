#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse

import numpy as np

dot = np.dot
norm = np.linalg.norm

mass_dict = {'H': 1.0, 'C': 12.0, 'O': 16.0, \
        'Pt': 195.0, 'Ag': 107.9}

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

def write_poscar(poscar, lattice, symbols, numbers, poses, fixes, vels):
    # write poscar
    content = 'POSCAR with Initialized Velocities\n'
    content += '{:>12.9f}\n'.format(1.0)
    for lat in lattice:
        content += '  {:>12.8f} {:>12.8f} {:>12.8f}\n'.format(*lat)
    content += ('  '+'{:>4s}'*len(symbols)+'\n').format(*symbols)
    content += ('  '+'{:>4d}'*len(numbers)+'\n').format(*numbers)
    content += 'Selective Dynamics\nDirect\n'
    for pos, fix in zip(poses,fixes):
        content += ('  '+'{:>16.12f}'*3+'{:>4s}'*3+'\n').format(*pos,*fix)
    content += '\n'
    for vel in vels:
        content += ('  '+'{:>16.12f}'*3+'\n').format(*vel)
                                                                           
    with open(poscar, 'w') as writer:
        writer.write(content)


def adjust_velocity(pa, pb, va, vb, mode):
    vec = pb - pa # a -> b
    vec = vec / norm(vec)

    # coef = 1.5 # enlarge
    if mode == 'a': # away
        va = -norm(va)*vec
        vb = norm(vb)*vec
    elif mode == 'c': # close
        va = norm(va)*vec
        vb = -norm(vb)*vec

    return va, vb


def initialize_velocity(poscar='POSCAR', aindices=[0,1], temperature=300):
    # parse internal coordinate
    coords = []
    for i in range(len(aindices)//3):
        coord = []
        coord.append(aindices[3*i])
        coord.append([aindices[3*i+1],aindices[3*i+2]])
        coords.append(coord)

    # first read poscar without velocity
    fname, scaling, lattice, symbols, numbers, poses, fixes = read_poscar('POSCAR')
    natoms = np.sum(numbers)

    elements = []
    for s, n in zip(symbols, numbers):
        elements.extend([s]*n)

    masses = []
    for e in elements:
        masses.append(mass_dict[e])
    masses = np.array(masses)

    nfrozens = 0
    for fix in fixes:
        for i in range(3):
            if fix[i] == 'T':
                nfrozens += 1

    # some units
    # UE, a.m. x A2 x fs2
    kb = 1.38064852e-23 # J / K
    JTOUE = 1/1.6605402*1e17 # J to a.m. x A2 / fs2
    BOLTOUE = 1/1.6605402*1.380648*1e-6 # a.m. x A2 / fs2 / K

    # generate velociteis using gaussian distribution
    vels = []
    for mass, fix in zip(masses, fixes):
        sigma = np.sqrt(BOLTOUE*temperature/mass)
        vel = np.zeros(3)
        for j in range(3):
            if fix[j] == 'T':
                vel[j] = np.random.normal(0,sigma,1)
        vels.append(vel)
    vels = np.array(vels)

    # give direction
    coef = 1.5
    for c in coords:
        mode = c[0]
        ia, ib = int(c[1][0])-1, int(c[1][1])-1
        va, vb = vels[ia], vels[ib]

        pa = np.dot(poses[ia], lattice)
        pb = np.dot(poses[ib], lattice)

        va, vb = adjust_velocity(pa, pb, va, vb, mode)
        vels[ia], vels[ib] = va, vb

    # remove drift
    for i in range(3):
        vels[:,i] = vels[:,i] - np.sum(vels[:,i]) / natoms

    # set fixes
    for i, fix in enumerate(fixes):
        for j in range(3):
            if fix[j] == 'F':
                vels[i][j] = 0.0

    # scale at temperature
    Ek = np.sum(1/2*dot(masses,vels**2))
    T = 2*Ek / (BOLTOUE*nfrozens) # current temperature
    scale = np.sqrt(temperature / T)
    vels = scale*vels

    Ek = np.sum(1/2*dot(masses,vels**2))
    T = 2*Ek / (BOLTOUE*nfrozens) # current temperature
    print('Current temperature is %d' %T)

    # writre poscar
    write_poscar(poscar+'-vel', lattice, symbols, numbers, poses, fixes, vels)
    print('Successfully write velocity included POSCAR.')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--poscar', nargs='?', default='POSCAR', \
            help='POSCAR File')
    parser.add_argument('-a', '--aindices', nargs='*', \
            required=True, help='Atom Indicies')
    parser.add_argument('-t', '--temperature', nargs='?', default=300, \
            type=float, help='Temperature')
    args = parser.parse_args()

    initialize_velocity(args.poscar, args.aindices, args.temperature)
