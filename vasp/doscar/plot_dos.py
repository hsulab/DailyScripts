#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse

import numpy as np

import matplotlib as mpl
mpl.use('Agg') #silent mode
from matplotlib import pyplot as plt

def plot_dos(dos_data_file, dos_fig_file):
    print(dos_data_file)
    with open(dos_data_file, 'r') as reader:
        lines = reader.readlines()
    lines = [line.strip().strip('#').split() for line in lines]

    col_names = lines[0]
    dos_data = np.array(lines[1:], dtype=float)

    energies = dos_data[:,0]

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8,8))
    ax.set_title('Projected Density of States', \
            fontsize=20, fontweight='bold')
    ax.set_xlabel('Energy / eV', fontsize=16)
    ax.set_ylabel('Density of States', fontsize=16)

    ax.plot(energies, len(energies)*[0], color='k')

    for col in range(1,dos_data.shape[1]):
        ax.plot(energies, dos_data[:,col], label=col_names[col])

    plt.legend()

    plt.savefig(dos_fig_file)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-f', '--datafile', \
            default='PDOS_USER.dat', help='PDOS Data File')
    parser.add_argument('-o', '--out', \
            default='DOS.png', help='DOS Figure')

    args = parser.parse_args()

    plot_dos(args.datafile, args.out)
