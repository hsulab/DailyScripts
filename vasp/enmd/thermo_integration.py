#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse

import numpy as np

from scipy import integrate

import argparse

import numpy as np

import matplotlib as mpl
mpl.use('Agg') #silent mode
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

def thermo_integration(datfile='BM-1.dat', nframe=-1):
    # read BM.dat file
    with open(datfile, 'r') as reader:
        lines = reader.readlines()
    lines = [line.strip().split() for line in lines if not line.startswith('#')]
    data = np.array(lines, dtype=float)

    print('Successfully read blue moon data from %s.' %datfile)

    coords = data[:nframe,1] # collective variable
    gradients = data[:nframe,5] # free energy gradients

    free_energies = []
    for i in range(len(gradients)):
        coord, grad = coords[:i+1], gradients[:i+1]
        v = integrate.trapz(grad, coord)
        free_energies.append(v)

    print('Successfully carry out the thermointegration.')

    # plot the figure
    fig, ax = plt.subplots(1, 1, figsize=(12,8))
    ax.set_title('Thermodynamic Integration', fontsize=20, fontweight='bold')

    a, = ax.plot(coords, free_energies, color='b', label='Free Energy')
    ax.set_xlabel('Collective Variable / a.u.', fontsize=16)
    # ax.set_xticks(np.arange(0, np.max(steps), 0.5))
    ax.xaxis.set_major_locator(MultipleLocator(0.1))
    ax.set_ylabel('Free Energy / eV', fontsize=16)

    ax2 = ax.twinx()
    b, = ax2.plot(coords, gradients, color='r', label='Free Energy Gradient')
    ax2.plot(coords, [0]*len(coords), c='k', ls='--')
    ax2.set_ylabel('Free Energy Gradient / a.u.', fontsize=16)

    plt.legend(handles=[a, b])

    plt.savefig('TI.png')

    print('Successfully plot the thermointegration process to TI.png.')

    # write free energies to file
    content = '# free energy\n'
    for step, force, energy in zip(coords, gradients, free_energies):
        content += '{:<12.4f}{:<12.4f}{:<12.4f}\n'.format(step, force, energy)

    with open('TI.dat', 'w') as writer:
        writer.write(content)

    print('Successfully write thermointegration data into TI.dat.')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file', required=True, \
            help='BM datfile')
    parser.add_argument('-nf', '--nframes', required=True, \
            type=int, help='Number of Frames')
    args = parser.parse_args()

    thermo_integration(args.file, args.nframes)
