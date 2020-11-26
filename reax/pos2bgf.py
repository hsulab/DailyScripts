#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse

import numpy as np

"""
Author: Jiayan Xu, jxu15@qub.ac.uk
Description:
    This script is for POSCAR/CONTCAR-BGF transformation.
Note:
    Atom connections are unnecessary for ReaxFF.
Updates:
    2019-12-13 LAMMPS use Y-along-B coordinate system.
"""


def read_poscar(poscar):
    # read poscar name
    poscar_name = os.path.basename(poscar)

    # only support VASP-5 format
    with open(poscar, 'r') as reader:
        lines = reader.readlines()
    lines = [line.strip().split() for line in lines]

    # read lattice constants
    scale = float(lines[1][0])
    abc = scale*np.array(lines[2:5], dtype=float)

    # read sys info
    elements = lines[5]
    numbers = lines[6]

    # read atoms
    atoms = []
    for element, number in zip(elements, numbers):
        atoms.extend([element]*int(number))
    natoms = len(atoms)

    # read coordinates
    frac_coords = np.array([line[:3] for line in lines[9:9+natoms]], dtype=float)
    cart_coords = np.dot(frac_coords, abc)

    return poscar_name, abc, atoms, frac_coords, cart_coords

norm = np.linalg.norm

def vec2deg(A, B):
    return np.arccos(np.dot(A,B)/(norm(A)*norm(B)))/np.pi*180

def pos2bgf(poscar):
    """Generate a lammps-like bgf file."""
    # read poscar
    poscar_name, abc, atoms, frac_coords, cart_coords = read_poscar(poscar)
    print('Read VASP5 %s' %os.path.basename(poscar))

    # FILTYP (A6,I5)
    content = 'BIOGRF 200\n' 

    # ('DESCRP',1X,A8)
    content += 'DESCRP ' + poscar_name + '\n' 

    # ('REMARK',1X,A)
    content += 'REMARK generated by pos2bgf.py\n' 

    # ('FORCEFIELD',1X,A8)
    # ('PERIOD',1X,3I1)

    # ('AXES',3X,A)
    content += 'AXES   xyz\n' 

    # ('CRYSTX',1X,6F11.5)
    A, B, C = abc[0], abc[1], abc[2]
    vol = np.dot(A, np.cross(B,C))
    a, b, c = norm(A), norm(B), norm(C)
    alpha, beta, gamma = vec2deg(B,C), vec2deg(A,C), vec2deg(A,B)

    # a along x
    lat_lmp = np.array([[a, 0, 0,], \
        [b*np.cos(gamma/180*np.pi), b*np.sin(gamma/180*np.pi), 0], \
        [0, 0, c]])

    a_lmp, b_lmp, c_lmp = lat_lmp[0], lat_lmp[1], lat_lmp[2]
    trans_matrix = 1/vol*np.dot(lat_lmp.T, \
            [np.cross(B,C),np.cross(C,A),np.cross(A,B)])

    lmp_coords = []
    for coord in cart_coords:
        lmp_coords.append(np.dot(trans_matrix, coord.T))
        #print(np.dot(trans_matrix, coord.T))
    lmp_coords = np.array(lmp_coords)

    content += ('CRYSTX ' + 6*'{:>11.5f}' + '\n')\
            .format(a, b, c, alpha, beta, gamma)
    cart_coords = lmp_coords
    
    # constraint
    #content += 'BOND RESTRAINT 1 2 0.50 7500.00 7500.00 0.0000000       0       0\n'
    
    # atoms
    # ('ATOM'|'HETATM',1X,I5,1X,A5,1X,A3,1X,A1,1X,A5,3F10.5,1X,A5,I3,I2,1X,F8.5)  
    atom_format = 'HETATM {:>5d} {:>5s} {:>3s} {:>1s} {:>5s}'+3*'{:>10.5f}'+\
            ' {:>5s}{:>3d}{:>2d} {:>8.5f}\n'
    for i, (atom, cart_coord) in enumerate(zip(atoms, cart_coords)):
        content += atom_format.format(i+1, atom, '', '', '', \
                cart_coord[0], cart_coord[1], cart_coord[2], atom, 0, 0, 0)

    # END
    content += 'END\n'

    bgf_path = os.path.abspath(poscar)+'.bgf'
    with open(bgf_path, 'w') as writer:
        writer.write(content)

    print('Write to BGF file %s.' %os.path.basename(bgf_path))

    return 

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file', required=True, help='POSCAR/CONTCAR')
    #parser.add_argument('-o', '--out', help='POSCAR/CONTCAR')
    args = parser.parse_args()

    pos2bgf(args.file)
