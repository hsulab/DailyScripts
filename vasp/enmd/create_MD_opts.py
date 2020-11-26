#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import shutil

import argparse

import numpy as np

from vaspy.iter import AniFile

from ase import Atom, Atoms
from ase.io import read, write
from ase.io.vasp import write_vasp

"""
Author: Jiayan Xu
Description:
    Create vasp MD opts.
"""

class Step(object):
    def __init__(self, niter, temperature, energy):
        self.niter = int(niter)
        self.temperature = round(float(temperature), 1)
        self.energy = round(float(energy), 4)

def select_structure(datfile='MD.dat', region=[8000,20000], interval=1000):
    # read MD data file
    print('Read %s ...' %datfile)
    with open(datfile, 'r') as reader:
        lines = reader.readlines()
    lines = [line.strip().strip('#').split() for line in lines]

    comments = lines[0]
    data = np.array(lines[1:], dtype=float) # time temperature energy

    timesteps, temperatures, energies = data[:,0], data[:,1], data[:,2]

    start, end = region[0], region[1]

    # assert len(region) == 2 and 
    # find relatively low energy structure
    stepmins, emins = [], []
    nsections = int((end-start)/interval)
    for sec in range(nsections):
        s, e = start + sec*interval, start + (sec+1)*interval
        if sec + 1 == nsections:
            e = end
        step_min = np.argmin(energies[s:e]) + s
        e_min = energies[step_min]
        print('Time Step %d -> Potential Energy %.4f' %(step_min+1, e_min))
        stepmins.append(step_min) # be careful, python index
        emins.append(e_min)
    print('Finished ...')

    return stepmins, emins

def create_vasp_dirs(datfile='MD.dat', region=[8000,20000], \
        interval=1000, optdir='md-opts'):
    # check dir
    if os.path.exists(optdir):
        optdir = os.path.abspath(optdir)
        print(optdir)
    else:
        raise ValueError('Please create a directory')

    # check template
    tempdir = os.path.join(optdir, 'template')
    if os.path.exists(tempdir):
        print('template directory -> %s' %tempdir)
        for vaspfile in ['INCAR', 'POTCAR', 'KPOINTS', 'vasp.script']:
            if not os.path.exists(os.path.join(tempdir,vaspfile)):
                raise ValueError('No %s in %s' %(vapfile,tempdir))
    else:
        raise ValueError('Please create a template')

    # get structures
    stepmins, emins = select_structure(datfile, region, interval) 
    # be careful, python index

    # check
    if os.path.exists('OUT.ANI'):
        pass
    else:
        raise ValueError('No OUT.ANI')

    # read and convert
    pos_dict= {}
    outani = AniFile('OUT.ANI')
    for n, ani in enumerate(outani):
        if n in stepmins:
            print('Read Step %d from OUT.ANI' %(n+1))
            pos_dict[n] = ani.data

    atoms = read('POSCAR', format='vasp')

    # create dirs
    for stepmin, emin in zip(stepmins, emins):
        dirname = str(stepmin+1).zfill(6) # outani index
        dirpath = os.path.abspath(os.path.join(optdir,dirname))
        print(dirpath)
        if os.path.exists(dirpath):
            ifremove = input('Remove %s [yes/no]: ' %dirpath)
            if ifremove == 'yes':
                shutil.rmtree(dirpath)
        # copy INCAR, KPOINTS, POTCAR, vasp.script
        shutil.copytree(tempdir, dirpath)
        print('--> Create Directory')
        atoms.set_positions(pos_dict[stepmin])
        # poscar
        poscar_path = os.path.join(dirpath, 'POSCAR')
        write(poscar_path, atoms, direct=True)
        with open(poscar_path, 'r') as reader:
            lines = reader.readlines()
        lines[0] = 'Step %d Energy %.4f\n' %(stepmin+1, emin)
        with open(poscar_path, 'w') as writer:
            writer.writelines(lines)
        print('--> Write POSCAR')
        # vasp.script
        vs_path = os.path.join(dirpath, 'vasp.script')
        with open(vs_path, 'r') as reader:
            lines = reader.readlines()
        lines[1] = '#PBS -N MD-%d\n' %(stepmin+1)
        with open(vs_path, 'w') as writer:
            writer.writelines(lines)

    return

if __name__ == '__main__':
    parser = argparse.ArgumentParser() 
    parser.add_argument('-m', '--mode', choices=('s', 'c'), \
            required=True, help='select or create')
    parser.add_argument('-f', '--file', help='MD.dat file')
    parser.add_argument('-r', '--region', required=True, \
            nargs='+', type=int, help='MD steps region')
    parser.add_argument('-i', '--interval', required=True, \
            type=int, help='selection interval')
    args = parser.parse_args()

    if args.mode == 's':
        select_structure(args.file, args.region, args.interval)
    elif args.mode == 'c':
        create_vasp_dirs(args.file, args.region, args.interval)

    pass
