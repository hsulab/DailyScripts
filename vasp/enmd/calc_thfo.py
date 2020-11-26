#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse

import numpy as np

dot = np.dot
norm = np.linalg.norm
inv = np.linalg.pinv # pseudo inverse

import matplotlib as mpl
mpl.use('Agg') #silent mode
from matplotlib import pyplot as plt

"""
Author: Jiayan Xu
"""

MASS_DICT = {'H': 1.0, 'C': 12.0, 'O': 16.0, \
        'Pt': 195.0, 'Ag': 107.9}

class ReactCoord():
    def __init__(self, tag, aindices, rctyp, coefs=None):
        self._tag = tag
        self.aindices = aindices
        self._rcval = 0
        self.rctyp = rctyp
        self.coefs = coefs

    @property
    def tag(self):
        return self._tag

    @property
    def rcval(self):
        return self._rcval

    @rcval.setter
    def rcval(self, rcval):
        self._rcval = rcval
        return
    
    def __repr__(self):
        curformat = 'TAG: {:>2s} ATOMS: ' + '{:<4d}'*len(self.aindices)
        content = curformat.format(self.tag, *self.aindices)
        return content

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


def read_thfo(filename='fort.129', natoms=2, nframes=3):
    # get file object
    reader = open(filename, 'r')
        
    nstep = 0
    frames = []
    while True:
        # find MD step
        line = reader.readline()
        if line.startswith('--------------------    MD STEP'):
            # count and reset data
            nstep += 1
            cur_frame = []
            irbm = False # If Read Blue Moon

        # read section
        if line.startswith(' -> CURRENT POSITIONS, '):
            for i in range(natoms):
                cur_atom = []
                line = reader.readline().strip().split()
                cur_atom.extend(line) # positions, velocities
                line = reader.readline().strip().split()
                cur_atom.extend(line) # forces
                cur_frame.append(cur_atom)
            frames.append(cur_frame)
            irbm = True

        # check step
        if nstep == nframes and irbm == True:
            break

        if not line:
            break
    reader.close()

    # check output
    frames = np.array(frames, dtype=float)
    poses, vels, forces = [], [], []
    for frame in frames:
        poses.append(frame[:,0:3])
        vels.append(frame[:,3:6])
        forces.append(frame[:,6:9])
    poses = np.array(poses) # direct coordinates
    vels = np.array(vels)
    forces = np.array(forces)

    return poses, vels, forces

def calc_thfo(thfofile, nframes, rcs, para, potim):
    # general setting
    nrc = len(rcs) # number of reactive coordinates

    # read poscar ...
    fname, scaling, lattice, \
            symbols, numbers, atoms, poses, fixes = read_poscar('POSCAR')
    natoms = np.sum(numbers)

    masses = []
    for atom in atoms:
        for i in range(3):
            masses.append(MASS_DICT[atom])
    masses = np.array(masses)
    masses_inv = 1.0 / masses

    UT, UL = 1e-15, 1e-10
    AMTOKG, EVTOJ = 1.6605402e-27, 1.60217733E-19
    FACT = AMTOKG*UL/UT**2/(EVTOJ/UL)

    # read fort.129
    poses, vels, forces = read_thfo(thfofile, natoms, nframes)

    # jacobian matrix
    fegs, rcoords = [], []
    for nframe in range(nframes):
        print('--------------------    MD STEP %8d--------------------\n' \
                %(nframe+1))
        # info
        pos, vel, force = poses[nframe], vels[nframe], forces[nframe]
        cur_pos, cur_vel, cur_for = pos, vel.ravel(), force.ravel()

        # !!!!!! scale
        for i, force in enumerate(cur_for):
            cur_for[i] = cur_for[i] * masses_inv[i]
        if nframe == 0:
            las_mj = np.zeros(nrc,3*natoms)
            las_vel = np.zeros(3*natoms)
            las_for = np.zeros(3*natoms)

        # calculation
        print('-> REACTIVE COORDINATES\n')
        jacob = np.zeros((nrc,3*natoms))
        for i, rcd in enumerate(rcs):
            if rcd.tag == 'R':
                jacob[i,:] = jdis(rcd, natoms, cur_pos, lattice)
                print('{:>20s}{:>20.16f}'.format('  DISTANCE',rcd.rcval))
                #if rcd.rctyp != 0:
                #    jacob[i,:] = np.zeros((1,3*natoms))
            elif rcd.tag == 'S':
                jacob[i,:] = np.zeros((1,3*natoms))
                rcd.rcval = 0
                for j in range(i):
                    rcd.rcval += rcs[j].rcval*rcd.coefs[j]
                    jacob[i,:] += jacob[j,:]*rcd.coefs[j]
                    jacob[j,:] = np.zeros((1,3*natoms))
                print('{:>20s}{:>20.16f}'.format('  COMBINITION',rcd.rcval))
        print('\n', end='')
        rcoords.append(rcs[para].rcval)

        # Mxi
        mxi = calc_mxi(jacob, masses_inv, natoms, nrc)
        mj = dot(mxi,jacob) # Mxi dot Jacobian
        cur_mj = mj

        #print('las_mj', np.sum(las_mj))
        #print('cur_mj', np.sum(cur_mj))
        #print('las_for', np.sum(las_for))
        #print('cur_for', np.sum(cur_for))

        # finally free energy gradients calculation
        print('-> FREE ENERGY GRADIENT (INSTANTANEOUS FORCE)\n')
        print('{:>12s}  {:>20s}  {:>20s}  {:>20s}  '\
                .format('        RC', '    FEG', '    FEG1', '    FEG2'))
        if nframe >= 1:
            # print(cur_for)
            feg1 = -FACT*dot(cur_mj-las_mj,cur_vel.T)/potim
            feg2 = -0.25*dot(cur_mj+las_mj,(cur_for+las_for).T)
            #xx = 0
            #for j in range(3*natoms):
            #    xx += (cur_mj+las_mj)[3,j] * (cur_for+las_for)[j]
            #    print((cur_mj+las_mj)[3,j], (cur_for+las_for)[j], xx)
            #print('-----')
            #print(xx*-0.25)
            #for cf, lf in zip(cur_for, las_for):
            #    print(cf, lf)
            feg = feg1 + feg2
            # print((cur_mj+las_mj).shape)
        else:
            feg1 = np.zeros(nrc)
            feg2 = np.zeros(nrc)
            feg = np.zeros(nrc)
        fegs.append(feg)

        for i in range(nrc):
            print('{:>8s}{:>4d}  {:>20.16f}  {:>20.16f}  {:>20.16f}  '\
                    .format('      RC', i+1, feg[i], feg1[i], feg2[i]))
        print('\n', end='')

        # update las
        las_vel, las_for, las_mj = cur_vel.copy(), cur_for.copy(), cur_mj.copy()

    fegs = np.array(fegs)
    content = '#step gradient\n'
    for i, (rcoord, feg) in enumerate(zip(rcoords, fegs[:,para])):
        content += '%4d %8.4f %8.4f\n' %(i+1, rcoord, feg)
    with open('ffff.dat', 'w') as writer:
        writer.write(content)

    return fegs

def calc_jacob():
    pass

def calc_mxi(jacob, masses_inv, natoms, nrc):
    jacob_trans = jacob.T
    mijt = np.zeros((3*natoms,nrc)) # inv of mass matrix dot transpose jacob
    for i in range(3*natoms):
        for j in range(nrc):
            mijt[i,j] = masses_inv[i] * jacob_trans[i,j]
    mxi_inv = dot(jacob, mijt)
    mxi = inv(mxi_inv)

    return mxi

def jdis(rcd, natoms, dirposes, latt):
    """jacobian distance"""
    jdis = np.zeros(3*natoms)

    ia, ib = rcd.aindices[0], rcd.aindices[1]

    pa = dirposes[ia]
    pb = dirposes[ib]

    vect = dot(pa-pb,latt)
    dis = norm(vect)
    rcd.rcval = dis
    
    vect = vect / dis # norm

    jdis[3*ia:3*ia+3] = vect
    jdis[3*ib:3*ib+3] = -vect

    return jdis


if __name__ == '__main__':
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--ncons', type=int, \
            required=True, help='Number of Constraints')
    parser.add_argument('-nf', '--nframes', nargs='?', \
            default=None, type=int, help='Number of Steps to Average')
    parser.add_argument('-p', '--parameter', type=int, \
            required=True, help='Constrained Parameter Index')

    args = parser.parse_args()
    '''

    #read(filename='fort.129', natoms=1, nframes=3)
    # R_CO
    #para = 0
    #rcs = [
    #        ReactCoord('R', [0,1], 0)]
    # R_OPt
    #para = 0
    #rcs = [
    #        ReactCoord('R', [1,14], 0)]
    # R_CPt
    para = 0
    rcs = [
            ReactCoord('R', [0,17], 0)]
    # R_COt
    #para = 0
    #rcs = [
    #        ReactCoord('R', [0,2], 0)]
    # R_OPtL
    #para = 0
    #rcs = [
    #        ReactCoord('R', [1,16], 0)]
    # R_OPtR
    #para = 0
    #rcs = [
    #        ReactCoord('R', [1,15], 0)]
    # -----
    #para=2
    #rcs = [
    #        ReactCoord('R', [0,1], 1),
    #        ReactCoord('R', [1,14], 1),
    #        ReactCoord('S', [], 0, [1,-1]),
    #        ] 
    calc_thfo('fort.129', 300, rcs, para, 0.5)
