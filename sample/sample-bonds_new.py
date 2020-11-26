#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import shutil
import subprocess
import argparse

import numpy as np

norm = np.linalg.norm
inv = np.linalg.inv

# read and write files
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


# set lengths
def set_bond_length_top(lattice,poses,ia,top,length):
    # calc
    poses[ia] = poses[top].copy()
    poses[ia][2] = poses[top][2] + length/lattice[2][2]

    return poses

def set_bond_length_brg(lattice,poses,ia,brg,length):
    # need cartesian
    dir_poses = poses.copy()
    car_poses = np.dot(dir_poses,lattice)
    pos_o = car_poses[ia]

    centre = (car_poses[brg[0]] + car_poses[brg[1]])/2.

    dis_z = (length**2-(2.806/2.0)**2)**0.5
    pos_o = centre.copy()
    pos_o[2] = pos_o[2] + dis_z
    dir_pos_o = np.dot(pos_o,inv(lattice))
    dir_poses[ia] = dir_pos_o

    # check distance
    poses = dir_poses.copy()

    return poses

def set_bond_length_fcc(lattice,poses,ia,fcc,length):
    # need cartesian
    dir_poses = poses.copy()
    car_poses = np.dot(dir_poses,lattice)
    pos_o = car_poses[ia]

    centre = (car_poses[fcc[0]] + car_poses[fcc[1]])/2.
    centre[1] = centre[1] + 3**0.5/6.0*2.806

    dis_z = (length**2-(2.806*1.0/3**0.5)**2)**0.5
    pos_o = centre.copy()
    pos_o[2] = pos_o[2] + dis_z
    dir_pos_o = np.dot(pos_o,inv(lattice))
    dir_poses[ia] = dir_pos_o

    # check distance
    poses = dir_poses

    return poses

def set_bond_length_trimer(lattice,poses,ia,ib,ic,angle,length):
    # need cartesian
    dir_poses = poses.copy()
    car_poses = np.dot(dir_poses,lattice)
    pos_a = car_poses[ia]

    proj_y = length*np.sin(angle/2./180.*np.pi)
    proj_z = length*np.cos(angle/2./180.*np.pi)

    pos_b = pos_a.copy()
    pos_b[1] = pos_b[1] + proj_y
    pos_b[2] = pos_b[2] + proj_z
    dir_pos_b = np.dot(pos_b,inv(lattice))
    # print(dir_pos_b)

    pos_c = pos_a.copy()
    pos_c[1] = pos_c[1] - proj_y
    pos_c[2] = pos_c[2] + proj_z
    dir_pos_c = np.dot(pos_c,inv(lattice))
    # print(dir_pos_c)

    dir_poses[ib] = dir_pos_b
    dir_poses[ic] = dir_pos_c

    # check distance
    poses = dir_poses

    return poses

def set_bond_length(geo_file,sample_type,length,bonds):
    # read (direct) poses and lattice
    if geo_file == 'POSCAR':
        fname,scaling,lattice,components,numbers,poses,fixes = read_poscar(geo_file)
    elif geo_file == 'geo_in.gen':
        natoms,coord_type,components,elements, poses, lattice = read_gen(geo_file)
    else:
        raise ValueError('Unsupported Package.')

    # in and out must be (direct) poses
    if sample_type == 'bulk':
        assert len(bonds) == 2
        ia, ib = bonds[0], bonds[1]
        vec = poses[ia] - poses[ib]
        dis = norm(np.dot(lattice.T,vec))
        zoom = length/dis
        lattice = lattice * zoom
    elif sample_type == 'top': # also for dimer
        assert len(bonds) == 2
        ia, top = bonds[0], bonds[1]
        poses = set_bond_length_top(lattice,poses,ia,top,length)
    elif sample_type == 'brg':
        assert len(bonds) == 3
        ia, brg = bonds[0], bonds[1:]
        poses = set_bond_length_brg(lattice,poses,ia,brg,length)
    elif sample_type == 'fcc': # also for hcp
        assert len(bonds) == 4
        ia, fcc = bonds[0], bonds[1:]
        poses = set_bond_length_fcc(lattice,poses,ia,fcc,length)
    elif sample_type == 'trimer':
        assert len(bonds) == 4
        ia, ib, ic, angle = bonds[0], bonds[1], bonds[2], bonds[3]
        bonds = (ia,ib,ic,angle,length)
        poses = set_bond_length_trimer(lattice,poses,*bonds)
        bonds = bonds[:3]
    else:
        raise ValueError('Unsupported Smaple Type.')

    # write
    if geo_file == 'POSCAR':
        vels = []
        write_poscar(geo_file,scaling,lattice,components,numbers,poses,fixes,vels)
        elements = []
        for c, n in zip(components,numbers):
            elements.extend([c]*n)
    elif geo_file == 'geo_in.gen':
        write_gen(geo_file,natoms,coord_type,components,elements,poses,lattice)
    else:
        raise ValueError('Unsupported Package.')

    # run
    tol = 0.001
    for i in bonds[1:]:
        vec = poses[ia] - poses[i]
        dis = norm(np.dot(lattice.T,vec))
        assert np.fabs(dis-length) < tol

    return elements


# task related
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

def read_outcar(outcar='OUTCAR',**kwargs):
    # how many steps to read
    nframes = 10000
    natoms = 10000
    if kwargs:
        if 'nframes' in kwargs.keys():
            if kwargs['nframes'] > 0:
                nframes = int(kwargs['nframes'])

    # read OUTCAR
    frames = []
    fopen = open(outcar, 'r')
    count, flag = 0, True
    while flag:
        line = fopen.readline()
        if line.strip().startswith('POSITION'):
            fopen.readline() # segment line ---...---
            poses, forces = [], []
            for n in range(natoms):
                line = fopen.readline()
                if not line.strip().startswith('-----'):
                    data = line.strip().split()
                    poses.append(data[:3]) # x y z
                    forces.append(data[3:]) # fx fy fz
                else:
                    break
            else:
                raise ValueError('More thane 10k atoms, huh?')
            poses = np.array(poses, dtype=float)
            forces = np.array(forces, dtype=float)
            frames.append((poses, forces))
            count += 1
        if count == nframes or (not line):
            flag = False
    fopen.close()

    return frames

def read_detailedout(detailed='detailed.out'):
    frames = []

    fopen = open(detailed, 'r')
    while True:
        line = fopen.readline()
        if line.strip().startswith('Total Forces'):
            forces = []
            for n in range(10000):
                line = fopen.readline()
                if line.strip():
                    data = line.strip().split()
                    forces.append(data[1:4]) # x y z
                else:
                    break
            else:
                raise ValueError('More thane 10k atoms, huh?')
            forces = np.array(forces,dtype=float)
            frames.append((forces))
        if not line:
            break
    
    return frames


def get_poses_and_forces(software):
    # both in cartesian and in eV/Å
    if software == 'vasp':
        frames = read_outcar(outcar='OUTCAR',nframes=1)
        assert len(frames[0][0]) == len(frames[0][1])
        natoms = len(frames[0][0])
        positions, forces = frames[0][0], frames[0][1]
    elif software == 'dftb':
        natoms, coord_type, symbols,elements,poses,lattice = read_gen('geo_in.gen')
        positions = poses.copy()
        positions = np.dot(positions,lattice) # cartesian
        frames = read_detailedout(detailed='detailed.out')
        forces = frames[0] / 0.194469064593167E-01 # har/bohr to eV/AA
        assert natoms == len(positions) == len(forces)
        #print(positions)
        #print(forces)
    elif software == 'test':
        natoms = 3
        positions = [(0,0,0)]*3
        forces = [(0,0,0)]*3
    else:
        raise ValueError('Unsupported software.')


    return natoms, positions, forces

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
    status, output = subprocess.getstatusoutput(\
        "echo '\n\nBond Length %s Calculation Forward' >> print-out" % cur_dis)
    status, output = subprocess.getstatusoutput(run_software)
    if status:
        raise ValueError('Something wrong in calculation...')

    # back up current geo file
    geo_bak = os.path.basename(geo_file) + '_' + '%d' %(int(cur_dis*1000))
    shutil.copy2(geo_file,geo_dir+'/'+geo_bak)

    # get calculated total energy
    en = get_total_energy(software)

    return en

def differentiate_bond_distances(dat,geo,sys,bond_lengths,fd_step=0.01,test_on=False):
    """
    dat - dat file
    geo - file name of the geometry
    geo_dir - directory backing up calculated geometries
    """

    if geo == 'POSCAR':
        geo_dir = 'POSCARs'
        software = 'vasp'
    elif geo == 'geo_in.gen':
        software = 'dftb'
        geo_dir = 'GENs'
    else:
        raise ValueError('Unsupported Geo File')

    if test_on:
        software = 'test'

    if os.path.exists(geo_dir):
        shutil.rmtree(geo_dir)
    os.mkdir(geo_dir)

    # prepare
    with open(dat,'w') as writer:
        writer.write('# Length Energy\n')

    subprocess.getstatusoutput(\
        "echo 'Bond Length Calculation !!!' > print-out")

    # parse input
    sample_type, bonds = sys[0], sys[1]

    # run
    for length in bond_lengths:
        # current length
        length = round(length,1)
        print('Start Bond Length %.2f' %length)

        # run vasp
        for_dis = length + fd_step
        set_bond_length(geo,sample_type,for_dis,bonds)
        for_en = run_task(software,geo_dir,geo,for_dis)

        # run vasp
        bak_dis = length - fd_step
        set_bond_length(geo,sample_type,bak_dis,bonds)
        bak_en = run_task(software,geo_dir,geo,bak_dis)

        # calculate bond force
        bond_force = - (for_en-bak_en) / (for_dis-bak_dis)

        # save data
        with open(dat,'a') as adder:
            adder.write( ('{:>8.4}  '*3+'{:>12.8f} '*3+'\n')\
                    .format(length,for_dis,bak_dis,for_en,bak_en,bond_force) )

    return

def sample_forces(xyz,geo,sys,lengths,test_on=False):
    # check running package
    if geo == 'POSCAR':
        bak_dir = 'POSCARs'
        software = 'vasp'
        out_file = 'OUTCAR'
    elif geo == 'geo_in.gen':
        bak_dir = 'GENs'
        software = 'dftb'
        out_file = 'detailed.out'
    else:
        raise ValueError('Unsupported Geo File')

    if os.path.exists(bak_dir):
        shutil.rmtree(bak_dir)
    os.mkdir(bak_dir)

    for f in [xyz,'print-out']:
        with open(f,'w') as writer:
            writer.write('')

    if test_on:
        software = 'test'

    sample, bonds = sys[0], sys[1]

    for dis in lengths:
        elements = set_bond_length(geo,sample,dis,bonds)
        en = run_task(software,bak_dir,geo,dis)

        # back up current geo file
        out_bak = os.path.basename(out_file) + '_' + '%d' %(int(dis*1000))
        shutil.copy2(out_file,bak_dir+'/'+out_bak)

        natoms, pos, force = get_poses_and_forces(software)

        # write xyz
        content = '%d\n' %natoms
        content += 'energy=%12.8f  distance=%12.8f\n' %(en, dis)
        for e, p, f in zip(elements,pos,force):
            content += '  %s  ' %e
            content += ('  {:>12.6f}  '*6+'\n').format(*p,*f)
        content += '\n'

        with open(xyz,'a') as adder:
            adder.write(content)

    return

if __name__ == '__main__':
    # few setting
    sys = ('trimer',(0,1,2,180))
    geo = 'POSCAR'

    # lengths
    fd_step = 0.001
    lengths = np.linspace(0.8,2.4,33)
    bond_lengths = list(lengths)
    bond_lengths.extend(list(lengths+fd_step))
    bond_lengths.extend(list(lengths-fd_step))
    bond_lengths = [round(b,3) for b in bond_lengths]
    # print(bond_lengths)

    #differentiate_bond_distances(dat,geo,sys,bond_lengths,fd_step=0.01,test_on=True)
    xyz = 'bonds.xyz'
    sample_forces(xyz,geo,sys,bond_lengths,test_on=True)
