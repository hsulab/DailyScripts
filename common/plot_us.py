#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse

import numpy as np

import matplotlib as mpl
mpl.use('Agg') #silent mode
from matplotlib import pyplot as plt

from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
    AutoMinorLocator)

# few globals
FIGSIZE = (16,12)

LAYOUT = dict(\
    left=0.05,right=0.95, \
    bottom=0.10, top=0.85, \
    wspace=0.00,hspace=0.00)

TFONT = {'size': 28, 'weight': 'bold'}
LFONT = {'size': 24}

LINE = dict(linewidth=4,marker='o',markerfacecolor='w',markersize=10)

TITLE = r"PMF of $\bf{CO\ Desorption\ on\ Top\ Site}$ by Umbrella Sampling"
XLABEL = r"Reaction Coordinate $R_{C-Pt}$ (Å)"
YLABEL = "Free Energy (eV)"


def plot_pmf(ax,coords,energies):
    # box
    spwidth = 2
    ax.spines['top'].set_linewidth(spwidth)
    ax.spines['bottom'].set_linewidth(spwidth)
    ax.spines['left'].set_linewidth(spwidth)
    ax.spines['right'].set_linewidth(spwidth)

    # labels
    ax.set_title(TITLE, TFONT)
    ax.set_xlabel(XLABEL, LFONT)
    ax.set_ylabel(YLABEL, LFONT)

    # curves
    ax.plot([-100,100], [0,0], \
            color='k',linestyle='--',linewidth=4)

    # ax.scatter(coords, energies, color='r', marker='o', s=100,zorder=10)
    a, = ax.plot(coords,energies,color='r',zorder=100,**LINE,label='Free Energy')

    #ax.errorbar(coords, energies, errors, marker='s', mfc='g', \
    #    mec='g', ecolor='g', label='Error')
    # ticks
    ax.set_xlim(1.5,4.0)
    ax.xaxis.set_major_locator(MultipleLocator(0.5))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%4.1f'))
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))

    ax.set_ylim(-0.5,1.5)
    ax.yaxis.set_major_locator(MultipleLocator(0.5))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%4.1f'))
    ax.yaxis.set_minor_locator(MultipleLocator(0.1))

    #ax.tick_params(which='both', width=2)
    ax.tick_params(which='major', length=10, width=2, labelsize=20)
    ax.tick_params(which='minor', direction='in', width=2, length=5)

    return a,


def plot_us(datfile='fe.dat',pic='us.png'):
    # read data
    with open(datfile, 'r') as reader:
        lines = reader.readlines()
    lines = [line.strip().split() for line in lines \
            if not line.strip().startswith('#')]
    data = np.array(lines, dtype='float')

    coords = data[:,0] * 10.0 # nm to AA
    energies = data[:,1] / 96.485 # kJ/mol to eV
    errors = data[:,2]

    # figure
    fig, axarr = plt.subplots(nrows=1, ncols=1, figsize=FIGSIZE)
    # plt.suptitle('')

    # first, pmf
    a, = plot_pmf(axarr,coords,energies)

    plt.legend(handles=[a,],fontsize=24)

    fig.tight_layout()

    plt.savefig(pic)

if __name__ == '__main__':
    '''
    parser = argparse.ArgumentParser()

    parser.add_argument('-f', '--datafile', \
            default='PDOS_USER.dat', help='PDOS Data File')
    parser.add_argument('-o', '--out', \
            default='DOS.png', help='DOS Figure')

    args = parser.parse_args()
    '''

    plot_us()
