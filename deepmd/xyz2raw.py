#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import json
import argparse

from collections import Counter

from ase.io import read, write

from tqdm import tqdm

import dpdata 

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-df', '--datafile', 
        default='data.xyz', help='time series files'
    )
    parser.add_argument(
        '-ej', '--enjson', 
        default='e0.json', help='json file with single atom energies'
    )

    args = parser.parse_args()

    # 
    substract_baseline = False
    if os.path.exists(args.enjson):
        substract_baseline = True
        with open(args.enjson, 'r') as fopen:
            e0_dict = json.load(fopen)

    # sanity check, dpdata only needs species, pos, Z, force, virial 
    # ase-extxyz is inconsistent with quip-xyz, especially the force 
    frames = read(args.datafile, ':')
    print('number of frames ', len(frames))

    atomic_properties = ['numbers', 'positions', 'forces']
    calc_props = ['energy', 'forces']

    for atoms in tqdm(frames):
        # remove extra properties in atoms
        cur_properties = list(atoms.arrays.keys())
        for prop in cur_properties:
            if prop not in atomic_properties:
                #atoms.arrays.pop(prop)
                del atoms.arrays[prop]
        # atoms info 
        # del atoms.info['feature_vector']
        # TODO: check if calculator exists 
        atoms.calc = None # ase copys xyz info to SinglePointCalculator?
        stored_forces = atoms.arrays.get('forces', None)
        if stored_forces is not None:
            atoms.arrays['force'] = stored_forces.copy()
            del atoms.arrays['forces']

        # calc
        #cur_calc_props = list(atoms.calc.results.keys())
        #for prop in cur_calc_props:
        #    if prop not in calc_props:
        #        del atoms.calc.results[prop]
        # move forces to force

        # check e0
        if substract_baseline:
            chemical_symbols = atoms.get_chemical_symbols()
            sym_dict = Counter(chemical_symbols)
            tot_e0 = 0.0
            for elem, num in sym_dict.items():
                tot_e0 += num*e0_dict[elem]
            atoms.info['energy'] -= tot_e0
            #print(tot_e0)
            #print(sym_dict)

    write('dp_raw.xyz', frames)

    # 
    xyz_multi_systems = dpdata.MultiSystems.from_file(
        file_name='./dp_raw.xyz', 
        fmt='quip/gap/xyz'
    )
    print(xyz_multi_systems)
    xyz_multi_systems.to_deepmd_raw('./raw_data/')
    pass
