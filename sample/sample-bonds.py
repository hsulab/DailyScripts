#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import shutil
import subprocess

import numpy as np

norm = np.linalg.norm


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

    # print('Successfully read POSCAR, taking it as the reference...')

    return fname, scaling, lattice, symbols, numbers, poses, fixes

def write_poscar(poscar, scaling, lattice, symbols, numbers, poses, fixes, vels):
    # write poscar
    content = 'Auto\n'
    content += '{:>12.9f}\n'.format(scaling)
    for lat in lattice:
        content += '  {:>12.8f} {:>12.8f} {:>12.8f}\n'.format(*lat)

    # symbols and numbers
    content += ('  '+'{:>4s}'*len(symbols)+'\n').format(*symbols)
    content += ('  '+'{:>4d}'*len(numbers)+'\n').format(*numbers)

    # positions
    content += 'Selective Dynamics\nDirect\n'
    for pos, fix in zip(poses,fixes):
        content += ('  '+'{:>16.12f}'*3+'{:>4s}'*3+'\n').format(*pos,*fix)
    content += '\n'

    if vels:
        for vel in vels:
            content += ('  '+'{:>16.12f}'*3+'\n').format(*vel)
                                                                           
    # write
    with open(poscar, 'w') as writer:
        writer.write(content)

    return

def set_bond_length(poscar,ia,ib,length):

    fname, scaling, lattice, symbols, numbers, \
        poses, fixes = read_poscar(poscar)
    poses[ia][2] = poses[ib][2] + length/lattice[2][2]

    vels = []
    write_poscar(poscar,scaling,lattice,symbols,numbers,poses,fixes,vels)

    vec = poses[ia] - poses[ib]
    dis = norm(np.dot(lattice.T,vec))

    assert np.fabs(dis-length) < 1E-4 

    return

def get_total_energy():
    stat, en = \
        subprocess.getstatusoutput("grep 'sigma->' OUTCAR | tail -n 1")
    en = float(en.strip().split()[-1])

    return en

if __name__ == '__main__':
    # few setting
    ia, ib = 0, 15
    dat = 'bonds.dat'
    fd_step = 0.01

    poscars = './POSCARs'

    # prepare
    with open(dat,'w') as writer:
        writer.write('# Length Energy\n')

    if os.path.exists(poscars):
        shutil.rmtree(poscars)
    os.mkdir(poscars)

    subprocess.getstatusoutput(\
        "echo 'Bond Length Calculation !!!' > print-out")

    # run
    bonds = np.linspace(1.4,2.8,15)
    for length in bonds:
        # current length
        length = round(length,1)

        # run vasp
        for_dis = length + fd_step
        set_bond_length('POSCAR',ia,ib,for_dis)
        subprocess.getstatusoutput(\
                "echo '\n\nBond Length %s Calculation Forward' >> print-out" % for_dis)
        subprocess.getstatusoutput('aprun -n 24 $EXEC 2>&1 >> print-out')
        subprocess.getstatusoutput(\
                "cp POSCAR %s/POSCAR_%s" %(poscars,for_dis))

        for_en = get_total_energy()

        # run vasp
        bak_dis = length - fd_step
        set_bond_length('POSCAR',ia,ib,bak_dis)
        subprocess.getstatusoutput(\
                "echo '\n\nBond Length %s Calculation Backward' >> print-out" % bak_dis)
        subprocess.getstatusoutput('aprun -n 24 $EXEC 2>&1 >> print-out')
        subprocess.getstatusoutput(\
                "cp POSCAR %s/POSCAR_%s" %(poscars,bak_dis))

        bak_en = get_total_energy()


        # calculate bond force
        bond_force = - (for_en-bak_en) / (for_dis-bak_dis)

        # save data
        with open(dat,'a') as adder:
            adder.write( ('{:>8.4}  '*3+'{:>12.8f} '*3+'\n')\
                    .format(length,for_dis,bak_dis,for_en,bak_en,bond_force) )

