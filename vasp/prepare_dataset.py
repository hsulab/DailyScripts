#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os

import ase.io
from ase.io.vasp import read_vasp_out

MAXFRAME = 100000 # maximum steps in outcar, this is a very large number for safe

MAXLINE = 1000

def read_outcar(outcar='OUTCAR',verbose=True,wdat=False,**kwargs):
    # how many steps to read
    nframes = MAXFRAME
    if kwargs:
        if 'nframes' in kwargs.keys():
            if kwargs['nframes'] > 0:
                nframes = int(kwargs['nframes'])

    # open OUTCAR
    fopen = open(outcar, 'r')

    # first check number of atoms
    readLineCounter = 0
    while True:
        readLineCounter += 1
        line = fopen.readline()
        if line.startswith(' Dimension of arrays:'):
            line = fopen.readline()
            line = fopen.readline()
            natoms = int(line.strip().split()[-1])
            break
        if readLineCounter > MAXLINE:
            raise ValueError('No NIONS after reading %d lines.' %MAXLINE)
    #print(natoms)

    # read energy, free energy, positions, forces, stresses(if have)
    readLineCounter = 0
    results = {}
    if_en, if_pos = False, False
    while True:
        # read a single line
        line = fopen.readline()
        # read positions and forces
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
        # read energy and free energy
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

        # check if converged
        if line:
            if if_en and if_pos:
                readLineCounter += 1
                if_en, if_pos = False, False
                if readLineCounter == nframes:
                    break
        else:
            break
            # raise ValueError('position and energy not ')

    fopen.close()
    
    if verbose:
        print('Successfully read OUTCAR, get positions and forces ...')

    return frames, energy, free_energy


def out2xyz(basicName, basicPath, smpDirList, readIndex):
    #bulkPath = '../bulks/sampling'
    #smpDirList = ['300K','MP05','MP09','MP15','MP20'] # sampling path list
    #for dirPath in os.listdir(bulkPath):
    #    print(dirPath)
    for smpDir in smpDirList:
        outcarPath = os.path.abspath(os.path.join(basicPath, smpDir+'/OUTCAR'))
        smpXyzName = '%s-%s.xyz' %(basicName, smpDir.lower())
        atom_frames = read_vasp_out(outcarPath,index=readIndex)
        if not isinstance(atom_frames, list):
            atom_frames = [atom_frames]
        nAtomFrames = len(atom_frames)
        print('total number of frames %d' %(nAtomFrames))
        for frameCount, atoms in enumerate(atom_frames):
            atoms.info['description'] = '%s - %d/%d' \
                    %(smpXyzName, frameCount+1, nAtomFrames)

        ase.io.write(smpXyzName, atom_frames)

        print('write %s based on %s' %(smpXyzName, outcarPath))

        #exit()
        

if __name__ == '__main__':
    ## bulks
    #bulkPath = '../bulks/sampling'
    #smpDirList = ['300K','MP05','MP09','MP15','MP20'] # sampling path list

    #out2xyz('md-bulk', bulkPath, smpDirList, '0:100000')

    # vacancy
    bulkPath = '../vacancy/sampling'
    smpDirList = ['300K','MP05','MP09','MP15','MP20'] # sampling path list

    out2xyz('md-vacancy', bulkPath, smpDirList, '0:100000')

    exit()

    ## clusters
    #clusterPath = '../clusters/sampling'
    #smpDirList = ['Pt13', 'Pt19', 'Pt24', 'Pt38'] # sampling path list

    #out2xyz('md-cluster', clusterPath, smpDirList, '0:1000')

    # surfaces
    surfacePath = '../surfaces/opt'
    #smpDirList = ['100', '110', '111', '210', '211', '221', \
    #        '310', '311', '320', '321', '322', '331', '332'] # sampling path list

    #for smpDir in smpDirList:
    #    vaspoutPath = os.path.abspath(os.path.join(surfacePath, smpDir+'/vasp.out'))
    #    with open(vaspoutPath, 'r') as reader:
    #        lines = reader.readlines()
    #    if lines[-1] != ' reached required accuracy - stopping structural energy minimisation\n':
    #        print(vaspoutPath)

    # remove 210, cannot be optimized
    smpDirList = ['100', '110', '111', '211', '221', \
            '310', '311', '320', '321', '322', '331', '332'] # sampling path list

    out2xyz('sp-surface', surfacePath, smpDirList, -1)
