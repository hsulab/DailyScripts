#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys

import argparse

import numpy as np

norm = np.linalg.norm
inv = np.linalg.inv
dot = np.dot
cross = np.cross
pi = np.pi
sin = np.sin
cos = np.cos
arccos = np.arccos

description='''
Jiayan Xu
'''


def read_poscar(poscar='POSCAR', carform='vasp5'):
    """read POSCAR"""
    with open(poscar, 'r') as reader:
        lines = reader.readlines()
    lines = [line.strip().split() for line in lines]

    fname = ' '.join(lines[0]) # file description
    scaling = float(lines[1][0])
    lattice = np.array(lines[2:5], dtype=float)

    if carform == 'vasp5':
        symbols = lines[5]
        numbers = [int(i) for i in lines[6]]
        natoms = np.sum(numbers)
        dyntype = ' '.join(lines[7]) # dynamic type
        coorsys = lines[8] # coordinate system
        pstart = 9
    elif carform == 'vasp4':
        symbols = fname.split()
        numbers = [int(i) for i in lines[5]]
        natoms = np.sum(numbers)
        dyntype = ' '.join(lines[6]) # dynamic type
        coorsys = lines[7] # coordinate system
        pstart = 8
    else:
        raise ValueError('Unknown POSCAR/CONTCAR format.')

    # read positions
    poses, fixes = [], []
    for coord in lines[pstart:pstart+natoms]:
        poses.append(coord[:3])
        fixes.append(coord[3:])
    poses = np.array(poses, dtype=float)

    print('Successfully read POSCAR in %s format.' %carform)

    return fname, scaling, lattice, symbols, numbers, poses, fixes


def adjust_poses(poses, refposes):
    """poses, poses_ref should be in direct"""
    # atoms, take last step as reference
    for i in range(len(poses)):
        for x in range(3):
            move = round(poses[i][x] - refposes[i][x], 0)
            poses[i][x] -= move
    refposes = poses.copy()

    return poses, refposes


def write_cif(poscar, cif_name):
    """"""
    fname, scaling, lattice, symbols, numbers, refposes, fixes \
            = read_poscar('POSCAR')

    # read poscar/contcar
    fname, scaling, lattice, symbols, numbers, poses, fixes = read_poscar(poscar)

    a, b, c = lattice
    lx, ly, lz = norm(a), norm(b), norm(c) # angstrom

    vol = dot(a, cross(b, c))

    alpha = arccos(dot(b,c)/ly/lz)/pi*180
    beta = arccos(dot(a,c)/lx/lz)/pi*180
    gamma = arccos(dot(a,b)/lx/ly)/pi*180 # degree

    # for non-orthogonal lattice, a along x or b along y is a problem
    lat_ax = np.array([\
            [lx, 0, 0,], \
            [ly*cos(gamma/180*pi), ly*sin(gamma/180*pi), 0], \
            [0, 0, lz]\
            ])

    lat_by = np.array([\
            [lx*sin(gamma/180*pi), lx*cos(gamma/180*pi), 0,],
            [0., ly, 0.],
            [0., 0., lz]\
            ])

    #a_ax, b_ax, c_ax = lat_ax[0], lat_ax[1], lat_ax[2]
    #trans_matrix = 1/vol*dot(lat_ax.T,[cross(b,c),cross(c,a),cross(a,b)])

    a_by, b_by, c_by = lat_by[0], lat_by[1], lat_by[2]
    trans = 1/vol*dot(lat_by.T,[cross(b,c),cross(c,a),cross(a,b)])

    poses, refposes = adjust_poses(poses, refposes)
    cartposes = dot(poses, lattice)
    cartposes = dot(cartposes, trans.T)
    poses = dot(cartposes, inv(lat_by))

    # write cif
    content = ''

    firstline = 'data_'
    for symbol in symbols:
        firstline += symbol
    firstline += '\n'
    content += firstline

    #
    content += '{:<30}{:<20}\n'.format('_audit_creation_method', '\'pos2cif.py\'')
    content += '{:<30}{:<20.16f}\n'.format('_cell_length_a', lx)
    content += '{:<30}{:<20.16f}\n'.format('_cell_length_b', ly)
    content += '{:<30}{:<20.16f}\n'.format('_cell_length_c', lz)
    content += '{:<30}{:<20.1f}\n'.format('_cell_angle_alpha', alpha)
    content += '{:<30}{:<20.1f}\n'.format('_cell_angle_beta' , beta)
    content += '{:<30}{:<20.1f}\n'.format('_cell_angle_gamma', gamma)
    content += '{:<30}{:<20}\n'.format('_symmetry_space_group_H-M' , '\'P1\'')
    content += '{:<30}{:<20}\n'.format('_symmetry_Int_Tables_number' , '\'1\'')
    content += '{:<30}{:<20}\n'.format('_symmetry_cell_setting' , '\'triclinic\'')
    content += '{:<30}\n'.format('loop_')
    content += '{:<30}\n'.format('_symmetry_equiv_pos_as_xyz')
    content += '{:<30}\n'.format('x,y,z')

    content += '\n{:<30}\n'.format('loop_')
    content += '{:<30}\n'.format('_atom_site_label')
    content += '{:<30}\n'.format('_atom_site_type_symbol')
    content += '{:<30}\n'.format('_atom_site_occupancy')
    content += '{:<30}\n'.format('_atom_site_fract_x')
    content += '{:<30}\n'.format('_atom_site_fract_y')
    content += '{:<30}\n'.format('_atom_site_fract_z')
    content += '{:<30}\n'.format('_atom_site_U_iso_or_equiv')

    count = 0
    for symbol, number in zip(symbols,numbers):
        for i in range(number):
            content += '{:<5s}{:<5s}{:<8s}{:<20.12f}{:<20.12f}{:<20.12f}{:<8s}\n'\
                .format(symbol+str(i+1),symbol,str(1.0000),\
                *poses[count],str(0.0000))
            count += 1

    with open(cif_name, 'w') as f:
        f.write(content)

    print('Successfully write %s to %s.' %(poscar, cif_name))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-v', '--version', action='version', \
        version='%(prog)s 1.1')
    parser.add_argument('-p', '--poscar', nargs='?', \
        default='POSCAR', help='POSCAR/CONTCAR')
    parser.add_argument('-cn', '--cifname', nargs='?', \
        default=None, help='Number of Frames')

    args = parser.parse_args()

    if not args.cifname:
        cifname = os.path.basename(os.getcwd()) + '.cif'

    write_cif(args.poscar, cifname)
