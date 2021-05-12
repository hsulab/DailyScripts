#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse

import numpy as np

from ase.io import read, write

import matplotlib as mpl
mpl.use('Agg') #silent mode
from matplotlib import pyplot as plt


def rms_dict(x_ref, x_pred):
    """ Takes two datasets of the same shape and returns a dictionary containing RMS error data"""

    x_ref = np.array(x_ref)
    x_pred = np.array(x_pred)

    if np.shape(x_pred) != np.shape(x_ref):
        raise ValueError('WARNING: not matching shapes in rms')

    error_2 = (x_ref - x_pred) ** 2

    average = np.sqrt(np.average(error_2))
    std_ = np.sqrt(np.var(error_2))

    return {'rmse': average, 'std': std_}

def energy_plot(ener_in, ener_out, ax, title='Plot of energy'):
    """ Plots the distribution of energy per atom on the output vs the input"""
    # check
    ener_in, ener_out = np.array(ener_in), np.array(ener_out)

    # scatter plot of the data
    ax.scatter(ener_in, ener_out)

    # get the appropriate limits for the plot
    for_limits = np.concatenate((ener_in,ener_out), axis=None)
    elim = (for_limits.min() - 0.05, for_limits.max() + 0.05)
    ax.set_xlim(elim)
    ax.set_ylim(elim)

    # add line of slope 1 for refrence
    ax.plot(elim, elim, c='k')

    # set labels
    ax.set_xlabel('Reference / eV')
    ax.set_ylabel('Prediction / eV')

    #set title
    ax.set_title(title)

    # add text about RMSE
    _rms = rms_dict(ener_in, ener_out)
    rmse_text = 'RMSE:\n' + str(np.round(_rms['rmse'], 3)) + ' +- ' + str(np.round(_rms['std'], 3)) + 'eV/atom'
    ax.text(0.9, 0.1, rmse_text, transform=ax.transAxes, fontsize='large', \
            horizontalalignment='right', verticalalignment='bottom')

    return 

def force_plot(in_force, out_force, ax, symbol='H', title='Plot of force'):
    """ Plots the distribution of force components per atom 
        on the output vs the input only plots for the given atom type(s)"""
    # extract data for only one species
    in_force = in_force[symbol]
    out_force = out_force[symbol]

    # scatter plot of the data
    ax.scatter(in_force, out_force)

    # get the appropriate limits for the plot
    for_limits = np.array(in_force + out_force)
    flim = (for_limits.min() - 1, for_limits.max() + 1)
    ax.set_xlim(flim)
    ax.set_ylim(flim)

    # add line of
    ax.plot(flim, flim, c='k')

    # set labels
    ax.set_ylabel('force by GAP / (eV/Å)')
    ax.set_xlabel('force by VASP / (eV/Å)')

    #set title
    ax.set_title(title)

    # add text about RMSE
    _rms = rms_dict(in_force, out_force)
    rmse_text = 'RMSE:\n' + str(np.round(_rms['rmse'], 3)) + ' +- ' + str(np.round(_rms['std'], 3)) + 'eV/Å'
    ax.text(0.9, 0.1, rmse_text, transform=ax.transAxes, fontsize='large', horizontalalignment='right',
        verticalalignment='bottom')


def calc_free_energy():
    pass


def extract_energy_and_forces(atom_frames,calc=None,atomic=True):
    """
    Electronic free energy and Hellman-Feynman forces
    """
    energies_dft, forces_dft = [], {}
    energies_gap, forces_gap = [], {}

    for atoms in atom_frames: # free energy per atom
        # basic info
        symbols = atoms.get_chemical_symbols()
        if atomic:
            natoms = len(atoms)
        else:
            natoms = 1
        # energy
        free_energy = atoms.get_potential_energy(force_consistent=True) # electronic free energy
        energies_dft.append(free_energy/natoms)
        # force
        forces = atoms.get_forces()
        for sym, force in zip(symbols,forces):
            if sym in forces_dft.keys():
                forces_dft[sym].extend(force)
            else:
                forces_dft[sym] = list(force)
        if calc:
            # use quip to calculate gap predicted energy
            atoms.set_calculator(calc)
            free_energy = atoms.get_potential_energy(force_consistent=True) # electronic free energy
            energies_gap.append(free_energy/natoms)
            # force
            forces = atoms.get_forces()
            for sym, force in zip(symbols,forces):
                if sym in forces_gap.keys():
                    forces_gap[sym].extend(force)
                else:
                    forces_gap[sym] = list(force)

    return forces_dft, forces_gap, energies_dft, energies_gap

def read_arrays(datafile):
    with open(datafile, 'r') as reader:
        lines = reader.readlines()

    lines = [line.strip().split() for line in lines if not line.startswith('#')]

    data = np.array(lines, dtype=float)

    return data

if __name__ == '__main__':
    # parser
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '-df', '--datafiles', nargs='*', default=['en.dat'], 
        help='data files containing properties'
    )

    args = parser.parse_args()

    # check files
    if len(args.datafiles) == 1:
        # read energy 
        data = read_arrays(args.datafiles[0])
    
        # plot
        fig, axarr = plt.subplots(
            nrows=1, ncols=1, figsize=(12,8)
        )
        plt.suptitle('References vs. Predictions')

        ax = axarr 
        energy_plot(data[:,0], data[:,1], ax)

    else:
        raise ValueError('Unsupported data files...')

    plt.savefig('rms.png')
