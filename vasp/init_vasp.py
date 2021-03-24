#!/usr/bin/env python3

import os
import time
import json
import logging
import argparse

import shutil 

import numpy as np 

from ase.io import read, write
from ase.io.vasp import write_vasp 
from ase.calculators.vasp import Vasp2
from ase.constraints import FixAtoms

vasp_params = {
    # INCAR 
    "nwrite": 2, 
    "istart": 0, 
    "lcharg": False, 
    "lwave": False, 
    "lorbit": 10,
    "npar": 4,
    "xc": "vdw-df",
    "encut": 400,
    "prec": "Normal",
    "ediff": 1E-5,
    "nelm": 180, 
    "nelmin": 6, 
    "ispin": 1,
    "ismear": 1,
    "sigma": 0.2,
    "algo": "Fast", 
    "lreal": "Auto", 
    "isym": 0, 
    "ediffg": -0.05,
    "nsw": 200,
    "ibrion": 2,
    "isif": 2,
    "potim": 0.02, 
    # KPOINTS 
    "gamma": True,
    "kpts": [2, 2, 2], 
} 


def create_copt(atoms, indices):
    """"""
    # If constrained get distance and create fort.188
    if indices:
        if len(atom_idxs) > 2:
            raise ValueError("More than two atoms end with '_c'")
        pt1, pt2 = atoms.positions[indices[0]], atoms.positions[indices[1]]

        # Use Ax=b convert to cartisan coordinate
        distance = np.linalg.norm(pt1-pt2)

        # Create fort.188
        ts_content = '1\n3\n6\n4\n0.04\n%-5d%-5d%f\n0\n' % \
            (indices[0]+1, indices[1]+1, distance)

        with open('fort.188', 'w') as f:
            f.write(ts_content)

        content = 'Constrain TS -->\n'
        """ 
        content += "     fort.188 has been created.\n"
        content += '     ' + '-'*20 + '\n'
        content += "     atom number: {:<5d}{:<5d}\n".format(atom_idxs[0]+1, atom_idxs[1]+1)
        content += "     atom name: {} {}\n".format(*atom_names)
        content += "     distance: {:f}\n".format(distance)
        content += '     ' + '-'*20 + '\n'
        """ 

        # Set IBRION = 1
        #incar.set('IBRION', 1)
        #content += "{} -> {}\n".format("IBRION", "1")
        content += '     Note: IBRION must be 1.\n'

        content += '\n'

        return content
    else:
        return ''

def create_slurm():

    pass 


def create_vasp_inputs(atoms, cons_indices, directory='vasp-test'): 
    """only for molecules"""
    # set constraints in z-position 
    cons = FixAtoms(indices=cons_indices)
    atoms.set_constraint(cons) 

    # pseudo 
    pp_path = "/mnt/scratch/chemistry-apps/dkb01416/vasp/PseudoPotential"
    if 'VASP_PP_PATH' in os.environ.keys():
        os.environ.pop('VASP_PP_PATH')
    os.environ['VASP_PP_PATH'] = pp_path

    # vdw 
    vdw_envname = 'ASE_VASP_VDW'
    vdw_path = "/mnt/scratch/chemistry-apps/dkb01416/vasp/pot"
    if vdw_envname in os.environ.keys():
        _ = os.environ.pop(vdw_envname)
    os.environ[vdw_envname] = vdw_path
    # print(os.environ[vdw_envname])

    # check kpoints 
    kpts = np.linalg.norm(atoms.cell, axis=1).tolist()
    kpts = [int(20./k)+1 for k in kpts] 
    vasp_params['kpts'] = kpts 

    content = 'KPOINTS -->\n'
    content += "     Set k-point -> \033[1;35m{} {} {}\033[0m\n".format(*kpts)
    content += '\n'

    print(content)

    # write('POSCAR', atoms, format='vasp')
    calc = Vasp2(
        command = "vasp_std", 
        directory = directory,  
        **vasp_params, 
    )
    calc.write_input(atoms)

    ase_sort = os.path.join(directory, 'ase-sort.dat')
    #shutil.rmtree(ase_sort)
    os.remove(ase_sort)

    # rewrite poscar 
    poscar = os.path.join(directory, 'POSCAR')
    write_vasp(poscar, atoms, direct=True, vasp5=True)

    from collections import Counter 
    symbols = atoms.get_chemical_symbols() 
    all_atoms = Counter(symbols) 
    fixed_symbols = [symbols[i] for i in cons_indices]
    fixed_atoms = Counter(fixed_symbols)

    natoms = len(atoms) 
    nfixed = natoms - np.sum(list(fixed_atoms.values()))

    content = '\nPOSCAR -->\n'
    content += '    \033[4m Element       Numbers      \033[0m\n'
    for sym in all_atoms.keys():
        content += '         %2s \033[1;33m%4d\033[0m \033[32m%4d\033[0m(T)\033[31m%4d\033[0m(F)\n'\
                %(sym, all_atoms[sym], all_atoms[sym]-fixed_atoms[sym], fixed_atoms[sym])
    content += '    \033[4m                            \033[0m\n'
    content += '      %2s \033[1;33m%4d\033[0m \033[32m%4d\033[0m(T)\033[31m%4d\033[0m(F)\n'\
                %('Total', natoms, natoms-nfixed, nfixed)
    content += '\n'

    print(content)

    return 


if __name__ == '__main__': 
    # Set argument parser.
    parser = argparse.ArgumentParser()

    # Add optional argument.
    parser.add_argument("-t", "--task", nargs='?', const='SC', help="set task type")
    parser.add_argument("-k", "--kpoints", help="set k-points")
    parser.add_argument("-q", "--queue", choices=('Gold', 'Test'), help="pbs queue type")
    parser.add_argument("-n", "--ncpu", help="cpu number in total")

    parser.add_argument("-f", "--file", help="structure file in any format")
    parser.add_argument(
        "-z", "--zpos", type=float, default=4., 
        help="fixed atoms at bottom layer"
    )
    parser.add_argument(
        # "-c", "--copt", action='store_true', 
        "-c", "--copt", type=int, nargs=2, 
        help="use constrained optimisation for transition state search"
    )

    args = parser.parse_args()

    # Add all possible arguments in INCAR file.
    atoms = read(args.file)

    # sort atoms by symbols and z-positions especially for supercells 
    if True: 
        numbers = atoms.numbers 
        zposes = atoms.positions[:,2].tolist()
        sorted_indices = np.lexsort((zposes,numbers))
        atoms = atoms[sorted_indices]

    indices = [atom.index for atom in atoms if atom.position[2] < args.zpos]
    create_vasp_inputs(atoms, indices, './')

    if args.copt:
        pass 

    pass 
