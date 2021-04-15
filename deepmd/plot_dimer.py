#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
plot dimer curve 
"""

import argparse

import numpy as np

import matplotlib
matplotlib.use('Agg') #silent mode
import matplotlib.pyplot as plt

from ase.io import read, write 

def parse_dimer_frames(xyzfile):
    frames = read(xyzfile, ':16') 
    dimer_symbols = frames[0].get_chemical_symbols()

    data = []
    for atoms in frames:
        assert len(atoms) == 2 
        energy = atoms.get_potential_energy() 
        dist = np.linalg.norm(atoms[0].position-atoms[1].position)
        data.append([dist,energy])
    data = np.array(data) 

    return dimer_symbols, data 

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-df', '--datafile', \
            default='bias', help='time series files')

    args = parser.parse_args()

    #with open(args.datafile, 'r') as reader:
    #    lines = reader.readlines()

    #lines = [line.strip().split() for line in lines if not line.startswith('#')]
    symbols, data = parse_dimer_frames('./evaluated.xyz')

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12,8))
    ax.set_title(
        '%s Dimer Curve' %('-'.join(symbols)), 
        fontsize=20, 
        fontweight='bold'
    )
    
    ax.set_xlabel('Distance [Å]', fontsize=16)
    ax.set_ylabel('Energyr [eV]', fontsize=16)

    plt.plot(data[:,0], data[:,1])

    plt.savefig('dimer.png')
