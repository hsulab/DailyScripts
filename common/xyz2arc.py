#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import time 

import argparse

import numpy as np

"""
Author: Jiayan Xu
Description:
    Transfer OUT.ANI to OUT.arc for Material Studio Visualization.
    OUT.ANI from (modified) VASP and geo_end.xyz from DFTB+ can be converted.
"""
    
pi = np.pi
norm = np.linalg.norm
dot = np.dot
cross = np.cross
arccos = np.arccos

def read_lattice(fname):
    """POSCAR or geo_in.gen"""
    # read POSCAR for lattice 
    #with open(poscar, 'r') as reader:
    #    lines = reader.readlines()

    #mat = np.array([line.strip().split() for line in lines[2:5]], \
    #        dtype=float)

    with open(fname, 'r') as reader:
        lines = reader.readlines()

    natoms = int(lines[0].strip().split()[0])

    lat_idx = 2 + natoms + 1 # natoms, elements, atoms, origin, lattice

    mat = np.array([line.strip().split() for line in lines[lat_idx:lat_idx+3]], \
            dtype=float)

    return mat

def mat2lat(mat):
    """read lattice matrix in POSCAR and transform to bravis lattice"""
    # lattice vectors
    a, b, c = mat

    # volume of the cell
    vol = dot(a, cross(b, c))

    # unit A
    lx, ly, lz = norm(a), norm(b), norm(c)

    # unit degree
    alpha = arccos(dot(b,c)/ly/lz)/pi*180
    beta = arccos(dot(a,c)/lx/lz)/pi*180
    gamma = arccos(dot(a,b)/lx/ly)/pi*180

    # check if a along x of b along y
    lat_lmp = np.array([[lx, 0, 0,], \
            [ly*np.cos(gamma/180*np.pi), ly*np.sin(gamma/180*np.pi), 0], \
            [0, 0, lz]])
    
    a_lmp, b_lmp, c_lmp = lat_lmp[0], lat_lmp[1], lat_lmp[2]
    trans_matrix = 1/vol*np.dot(lat_lmp.T, \
            [np.cross(b,c),np.cross(c,a),np.cross(a,b)])

    return lx, ly, lz, alpha, beta, gamma, trans_matrix


def ani2arc(geo='POSCAR', ani='OUT.ANI', frame=None):
    # select frame
    if frame:
        frame = int(frame)

    # read lattice
    mat = read_lattice(geo)
    lx, ly, lz, alpha, beta, gamma, trans_matrix = mat2lat(mat)

    # number of atoms
    with open(ani, 'r') as reader:
        natoms = int(reader.readline())

    # read each frame
    count = 0 # number of frames
    content = '!BIOSYM archive 3\nPBC=ON\n'

    reader = open(ani, 'r') 
    for i, line in enumerate(reader):
        # read coordinates
        if i % (natoms+2) == 0:
            # write time
            cur_content = ('%80.4f\n' % 0.0)
            cur_content += '!DATE     %s\n' \
                    %time.asctime( time.localtime(time.time()))
            # write crystal
            cur_content += 'PBC%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f\n' %\
                    (lx, ly, lz, alpha, beta, gamma)
        elif 1 < i%(natoms+2):
            # write coordinates
            line = line.strip().split()
            atom = line[0]
            coord = np.array(line[1:4], dtype=float) # a long x or b along y
            coord = dot(trans_matrix, coord.T)
            cur_content += '%2s%16.9f%16.9f%16.9f%5s%2d%8s%8s%7.3f\n' %\
                    (atom, coord[0], coord[1], coord[2],
                    'XXXX', 1, 'xx', atom, 0.0)
        # add frame to content
        if (i+1) % (natoms+2) == 0:
            cur_content += 'end\nend\n'
            count += 1
            if frame and count == frame:
                content = '!BIOSYM archive 3\nPBC=ON\n' + cur_content
                print('Select Frame %d.' %count)
                break
            content += cur_content
            print('Write Frame %8d ...' %count, end='\r')
    print('Write Frame %8d ...' %count, end='\r')

    reader.close()

    # write arc to file
    if frame:
        arcname = os.path.basename(os.getcwd()) + '-%d.arc' %frame
    else:
        arcname = os.path.basename(os.getcwd()) + '.arc'
    with open(arcname, 'w') as f:
        f.write(content)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-gf", "--geoFile", default='geo_in.gen', \
            help="VASP/DFTB+ Geometry File")
    parser.add_argument("-xf", "--xyzFile", default='OUT.ANI', \
            help="VASP ANI FORMAT FILE or DFTB+ XYZ FORMAT FILE")
    parser.add_argument("-s", "--select", help="select given frame")

    args = parser.parse_args()
    ani2arc(args.geoFile, args.xyzFile, args.select)
