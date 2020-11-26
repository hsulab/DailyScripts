#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import shutil
import subprocess

import numpy as np

norm = np.linalg.norm
inv = np.linalg.inv

def read_gen(gen):
    with open(gen,'r') as reader:
        lines = reader.readlines()

    # natoms and coord type
    data = lines[0].strip().split()
    natoms = int(data[0])
    coord_type = data[1]
    if coord_type != 'F':
        raise ValueError('Only can F-coord be read.')

    # symbols
    symbols = lines[1].strip().split()

    # coord
    elements = []
    poses = []
    for line in lines[2:2+natoms]:
        data = line.strip().split()
        elements.append(symbols[int(data[1])-1])
        poses.append(data[2:])
    poses = np.array(poses,dtype=float)

    # lattice
    data = [line.strip().split() for line in lines[2+natoms:2+natoms+4]]
    data =  np.array(data,dtype=float)
    lat_org = data[0]
    lattice = data[1:]

    return natoms, coord_type, symbols, elements, poses, lattice

def write_gen(gen,natoms,coord_type,symbols,elements,poses,lattice):
    #
    sym_dict = {}
    for i, sym in enumerate(symbols):
        sym_dict[sym] = i+1
    #
    content = ('  {:<8d}  {:<s}\n').format(natoms,coord_type)
    content += ('{:>4s}'*len(symbols)+'\n').format(*symbols)
    for i, (sym, coord) in enumerate(zip(elements,poses)):
        content += ('{:>8d}  {:>4d}  '+'{:>16.12f}  '*3+'\n')\
                .format(i+1,sym_dict[sym],*coord)
    content += ('{:>16.12f}  '*3+'\n').format(*[0.]*3)
    for lat in lattice:
        content += ('{:>16.12f}  '*3+'\n').format(*lat)

    with open(gen,'w') as writer:
        writer.write(content)

    return


def set_bond_length_fcc(gen,ia,fcc,length):
    #
    natoms, coord_type, symbols, elements, dir_poses, lattice = read_gen(gen)

    # need cartesian
    car_poses = np.dot(dir_poses,lattice)
    pos_o = car_poses[ia]

    centre = (car_poses[fcc[0]] + car_poses[fcc[1]])/2.
    centre[1] = centre[1] + 3**0.5/6.0*2.806

    dis_z = (length**2-(2.806*1.0/3**0.5)**2)**0.5
    pos_o = centre.copy()
    pos_o[2] = pos_o[2] + dis_z
    dir_pos_o = np.dot(pos_o,inv(lattice))
    dir_poses[ia] = dir_pos_o

    # write poscar
    write_gen(gen,natoms,coord_type,symbols,elements,dir_poses,lattice)

    # check distance
    poses = dir_poses

    tol = 0.001
    vec = poses[ia] - poses[fcc[0]]
    dis = norm(np.dot(lattice.T,vec))
    assert np.fabs(dis-length) < tol

    vec = poses[ia] - poses[fcc[1]]
    dis = norm(np.dot(lattice.T,vec))
    assert np.fabs(dis-length) < tol

    vec = poses[ia] - poses[fcc[2]]
    dis = norm(np.dot(lattice.T,vec))
    assert np.fabs(dis-length) < tol

    return

def set_bond_length_top(gen,ia,top,length):
    #
    natoms, coord_type, symbols, elements, poses, lattice = read_gen(gen)

    poses[ia] = poses[top].copy()
    poses[ia][2] = poses[top][2] + length/lattice[2][2]

    # write poscar
    write_gen(gen,natoms,coord_type,symbols,elements,poses,lattice)

    vec = poses[ia] - poses[top]
    dis = norm(np.dot(lattice.T,vec))

    assert np.fabs(dis-length) < 1E-4 

    return


def get_total_energy(software):
    # check package
    if software == 'vasp':
        read_energy = "grep 'sigma->' OUTCAR | tail -n 1 | awk '{print $7}'"
    elif software == 'dftb':
        read_energy = "grep 'Extrapolated to 0' detailed.out | awk '{print $6}'"
    elif software == 'test':
        read_energy = "echo 100"
    else:
        raise ValueError('Unsupported software.')

    # return energy
    stat, en = subprocess.getstatusoutput(read_energy)
    en = float(en)

    return en

def run_task(software,geo_dir,geo_file,cur_dis):
    # check package
    if software == 'vasp':
        run_software = 'aprun -n 24 $EXEC 2>&1 >> print-out'
    elif software == 'dftb':
        run_software = 'aprun -cc none -n 2 -N 2 -S 1 -d 12 $EXEC 2>&1 >> print-out'
    elif software == 'test':
        run_software = 'echo test >> print-out'
    else:
        raise ValueError('Unsupported software.')

    # run!
    subprocess.getstatusoutput(\
        "echo '\n\nBond Length %s Calculation Forward' >> print-out" % cur_dis)
    subprocess.getstatusoutput(run_software)

    # back up current geo file
    geo_bak = os.path.basename(geo_file) + '_' + '%s' %(str(cur_dis))
    shutil.copy2(geo_file,geo_dir+'/'+geo_bak)

    return

if __name__ == '__main__':
    # few setting
    ia = 0
    top = 1

    dat = 'bonds.dat'
    fd_step = 0.01

    software = 'test'
    geo_file = 'geo_in.gen'
    geo_dir = './GENs'

    # prepare
    with open(dat,'w') as writer:
        writer.write('# Length Energy\n')

    if os.path.exists(geo_dir):
        shutil.rmtree(geo_dir)
    os.mkdir(geo_dir)

    subprocess.getstatusoutput(\
        "echo 'Bond Length Calculation !!!' > print-out")

    # run
    bonds = np.linspace(0.8,2.8,21)
    for length in bonds:
        # current length
        length = round(length,1)
        print('Start Bond Length %.2f' %length)

        # run vasp
        for_dis = length + fd_step
        set_bond_length_top(geo_file,ia,top,for_dis)
        run_task(software,geo_dir,geo_file,for_dis)
        for_en = get_total_energy(software)

        # run vasp
        bak_dis = length - fd_step
        set_bond_length_top(geo_file,ia,top,bak_dis)
        run_task(software,geo_dir,geo_file,bak_dis)
        bak_en = get_total_energy(software)

        # calculate bond force
        bond_force = - (for_en-bak_en) / (for_dis-bak_dis)

        # save data
        with open(dat,'a') as adder:
            adder.write( ('{:>8.4}  '*3+'{:>12.8f} '*3+'\n')\
                    .format(length,for_dis,bak_dis,for_en,bak_en,bond_force) )

