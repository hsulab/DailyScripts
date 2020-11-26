#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os
import sys

import argparse

import numpy as np

from ase import Atoms
from ase.build import make_supercell
from ase.io import read, write

help_context = """
Author: 
    jyxu
Description:
             __________ 
           /|          |
            |          |
       10 Å |  Vacuum  |
            |__________|
           \|____He____|
           /|          |
            |          |
       10 Å |  water   |
            |          |
           \|__________|
            |          |
            |          |
            |  bulk    |
            |          |
            |__________|
    
    1. Add water layer in selected region. (default 10 Å thick)
    2. Add Hellium slab to constrain water. (1 Å thick slab)
Notes:
    Water @ 298.15 K. 
    Rough estimation of volume.
    Add fixed Hellium (cub or hex) to constrain water.
"""


# ===== ====== ===== ===== ====== ===== ===== ====== ===== ===== 
# Water Molecules Number Calculation @ 298.15K
# ===== ====== ===== ===== ====== ===== ===== ====== ===== ===== 
def calc_water_num(volume):
    """Calculate number of water molecules at 298.15K with given volume"""
    # water density 
    water_molecule_weight = 18.0152 # g/mol
    water_density = 0.997074 # g/cm^3 298.15K
    n_avogadro = 6.02e23 # mol
    # volume = 16*13*10 # Å^3

    n_water_molecule_per_A = (water_density / water_molecule_weight) * n_avogadro * 1e-24
    
    return np.floor(n_water_molecule_per_A * volume)


# ===== ====== ===== ===== ====== ===== ===== ====== ===== ===== 
# Cubic or Hexagonal Box Generation
# ===== ====== ===== ===== ====== ===== ===== ====== ===== ===== 
def solve_plane_normal_vector(vec):
    """Set x=1, z=0"""
    y = -vec[0]/vec[1]
    
    return np.array([1., y, 0.])


def generate_constraints(boxtype, a, b, zrange):
    """generate a hexagonal / cubic box"""
    if boxtype == 'hex':
        constraints = '  # generate a hex box\n'
        # x direction
        constraints += '  # x-axis\n'
        vec = solve_plane_normal_vector(a) # plane normal vector
        d = np.dot(b, vec) / np.linalg.norm(vec) * np.linalg.norm(vec)
        constraints += '  over plane %.1f %.1f %.1f 0.\n' \
                %(vec[0], vec[1], vec[2])
        constraints += '  below plane %.1f %.1f %.1f %.1f\n' \
                %(vec[0], vec[1], vec[2], d)

        # y direction
        constraints += '  # y-axis\n'
        vec = solve_plane_normal_vector(b) # plane normal vector
        d = np.dot(a, np.array([1., 0., 0.]))
        constraints += '  over plane %.1f %.1f %.1f 0.\n' \
                %(vec[0], vec[1], vec[2])
        constraints += '  below plane %.1f %.1f %.1f %.1f\n' \
                %(vec[0], vec[1], vec[2], d)

        # z direction
        constraints += '  # z-axis\n'
        constraints += '  over plane 0. 0. 1. %.1f\n' %zrange[0] # bottom
        constraints += '  below plane 0. 0. 1. %.1f\n' %zrange[1] # top

        constraints += '  constrain_rotation x 0. 0. \n'
        constraints += '  constrain_rotation y 0. 0. \n'
        constraints += '  constrain_rotation z 0. 0. \n'
    elif boxtype == 'cub':
        constraints = '  # generate a cub box\n'
        # x direction 
        constraints += '  inside box 0. 0. %.1f %.1f %.1f %.1f\n' \
                %(zrange[0], np.linalg.norm(a), np.linalg.norm(b), zrange[1])

        constraints += '  constrain_rotation x 0. 0. \n'
        constraints += '  constrain_rotation y 0. 0. \n'
        constraints += '  constrain_rotation z 0. 0. \n'
    else:
        raise ValueError('Unknown boxtype for constraints generation.')

    return constraints 


# ===== ====== ===== ===== ====== ===== ===== ====== ===== ===== 
# Write structure part in a given format.
# ===== ====== ===== ===== ====== ===== ===== ====== ===== ===== 
def write_structure(strname, number, constraints):
    # start of structure
    structure = 'structure %s\n' %strname

    # number of molecules
    structure += '  number %d\n' %number

    # constraints 
    structure += constraints

    # end of structure
    structure += 'end structure\n'

    return structure 


# ===== ====== ===== ===== ====== ===== ===== ====== ===== ===== 
# ===== ====== ===== ===== ====== ===== ===== ====== ===== ===== 
def estimate_cluster_volume(poses, method='atomic'):
    if method == 'atomic':
        cluster_volume = len(poses)*4./3*np.pi*1.3**3 
    elif method == 'sphere':
        # sphere
        centre = np.zeros(3)
        for pos in poses:
            centre += pos / len(poses)
        radius = 0.
        for pos in poses:
            distance = np.linalg.norm(pos - centre)
            if distance > radius:
                radius = distance
        cluster_volume = 4./3.*np.pi*(radius+0.3)**3 # sphere
    elif method == 'MCE':
        # under development
        bonding_distance = 2.6 # Å
        connectivities = []
        for i in range(len(poses)):
            posi, connects = poses[i], []
            for j in range(1, len(poses)):
                posj = poses[j]
                if np.linalg.norm(posi-posj) <= bonding_distance:
                    connects.append(j)
            connectivities.append(connects)
            #print(connects)
    else:
        raise ValueError('Unknown method for cluster volume estimation.')

    return cluster_volume


def write_packmol_input(substrate, water, packmolin='autopackmol.inp'):
    """"""
    # general settings
    tolerance = 2.0 # Å
    filetype = 'xyz'
    water_thickness = 10.0 # Å
    cluster_metal = []

    # check file extension / format
    basename = os.path.basename(substrate).split('.')
    
    substrate_name, substrate_ext =  basename[0], basename[1]

    basename = os.path.basename(water).split('.')
    water_name, water_ext =  basename[0], basename[1]

    if water_ext != filetype:
        raise ValueError('Get me a %s-format water.' %filetype)

    try:
        atoms = read(substrate, format=substrate_ext)
    except:
        raise ValueError('substrate should either be in xyz format.')

    # adjust atom positions in z
    poses = atoms.get_scaled_positions()
    for n, pos in enumerate(poses):
        if pos[2] > 0.90:
            poses[n][2] = pos[2]-1.0
    print(poses)
    atoms.set_scaled_positions(poses)

    if substrate_ext != filetype:
        substrate = substrate_name + '.' + filetype
        write(substrate, atoms, format=filetype)

    # lattice type detection
    a, b, c = atoms.get_cell()
    if np.fabs(np.dot(a, c)) < 1e-6 and np.fabs(np.dot(b, c)) < 1e-6:
        if np.fabs(np.dot(a, b)) < 1e-6:
            boxtype = 'cub'
        else:
            boxtype = 'hex'
    else:
        raise ValueError('Only support lattice Z perpendicular to XY plane.')

    # write packmol.inp
    content = '# Created by generate_packmol_in.py\n'

    # molecule tolerance
    content += '# molecules will be separated at least %s Å\n' %tolerance
    content += 'tolerance %.1f\n' %tolerance
    content += '\n'

    # file type
    content += 'filetype %s\n' %filetype
    content += 'output %s\n' %(substrate_name + '_' + water_name \
              + '.' + filetype)
    content += '\n'

    # add structures
    structures = '# add %f Å thick water layers on the substrate\n' %water_thickness

    # substrate
    posz = np.ceil(np.max([z for [x, y, z], s in \
            zip(atoms.get_positions(), atoms.get_chemical_symbols())]))

    constraints = generate_constraints(boxtype, a, b, [0.,posz])
    structures += write_structure(substrate, 1, constraints)
    structures += '\n'

    # water layer, add shield for cluster if neccessary
    cluster_poses = [pos for pos, s in \
            zip(atoms.get_positions(), atoms.get_chemical_symbols()) \
            if s in cluster_metal]
    if len(cluster_poses) != 0:
        volume = - estimate_cluster_volume(cluster_poses, 'atomic')
        posz_sub = np.floor(np.min([z for [x, y, z], s in \
            zip(atoms.get_positions(), atoms.get_chemical_symbols())\
            if s in cluster_metal]))
    else:
        volume = 0.
        posz_sub = posz

    constraints = generate_constraints(boxtype, a, b, [posz_sub,posz_sub+water_thickness])

    volume += np.dot(np.cross(a, b), [0., 0., water_thickness])
    water_num = calc_water_num(volume)

    structures += write_structure(water, water_num, constraints)
    structures += '\n'

    # add hellium slab
    if boxtype == 'cub':
        He_cub = Atoms('He2', scaled_positions=[[0.,0.,0.], [0.5, 0.5, 0.]], \
                cell=[[4.242, 0., 0.], [0., 4.242, 0.], [0., 0., 1.0]])
        superx, supery = round(np.linalg.norm(a) / 4.242, 0), \
                round(np.linalg.norm(b) / 4.242, 0)
        He_slab = make_supercell(He_cub, [[superx,0,0],[0,supery,0],[0,0,1]])
    elif boxtype == 'hex':
        He_hex = Atoms('He', scaled_positions=[[0.3333, 0.6667, 0.]], \
                cell=[[3.161, -1.825, 0.], [0., 3.650, 0.], [0., 0., 1.0]])
        superx, supery = np.linalg.norm(a) / 3.650, np.linalg.norm(b) / 3.650
        He_slab = make_supercell(He_hex, [[superx,0,0],[0,supery,0],[0,0,1]])
    write('He_slab.xyz', He_slab, format='xyz')
    constraints = generate_constraints(boxtype, a, b, \
            [posz_sub+water_thickness, posz_sub+water_thickness+1.0])
    structures += write_structure('He_slab.xyz', 1, constraints)
    structures += '\n'

    content += structures
    
    # write inp 
    print(content)
    with open(packmolin, 'w') as writer:
        writer.write(content)
    print('Successfully generate %s for packmol input.' %packmolin)

    # under development, autorun packmol (which is not security)
    #packmol = '../packmol/packmol'
    #os.popen(packmol)


if __name__ == '__main__':
    parser = argparse.ArgumentParser() 
    parser.add_argument('-s', '--substrate', default='substrate.xyz', \
            help='substrate')
    parser.add_argument('-w', '--water', default='water.xyz', \
            help='water')

    args = parser.parse_args()

    write_packmol_input(args.substrate, args.water)
