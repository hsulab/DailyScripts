#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
plot time versus colvar
"""

import argparse
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use('Agg') #silent mode
import matplotlib.pyplot as plt

from ase.io import read, write

def plot_status(ax, data: list):
    """"""
    # unpack data
    energies, maxforces, varforces = data

    nsteps = len(energies) # number of optimisation steps
    assert nsteps == maxforces.shape[0] == varforces.shape[0]

    # forces
    ax.set_xlabel('Number of Step', fontsize=16)
    ax.set_ylabel('Force [eV/AA]', fontsize=16)

    ax.plot(range(nsteps), maxforces)
    ax.fill_between(
        range(nsteps), 
        maxforces-np.sqrt(varforces), maxforces+np.sqrt(varforces), 
        alpha=0.2
    )

    # energies
    ax2 = ax.twinx()
    ax2.set_ylabel('Potential Energy [eV]', fontsize=16)
    ax2.plot(range(nsteps), energies)

    return

def read_outcar(outcar):
    """read outcar and return energies and forces"""
    frames = read(outcar, ':')
    energies = [atoms.get_potential_energy() for atoms in frames]
    force_arrays = np.array([atoms.get_forces().flatten() for atoms in frames])
    print(force_arrays.shape)
    maxforces = np.max(force_arrays, axis=1)
    print(maxforces.shape)
    varforces = np.var(force_arrays, axis=1)
    print(varforces.shape)

    print("maxforce: ", maxforces)

    return energies, maxforces, varforces

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-d', '--dirs', nargs='*', 
        default=['./'], help='vasp directories'
    )
    args = parser.parse_args()

    outcars = [Path(d) / 'OUTCAR' for d in args.dirs]
    print(outcars)

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12,8))
    ax.set_title(
        'Optimisation Status', 
        fontsize=24, fontweight='bold'
    )

    for outcar in outcars:
        data = read_outcar(outcar)
        plot_status(ax, data)

    plt.savefig('fopt.png')
