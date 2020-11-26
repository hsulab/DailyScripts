#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse

import numpy as np

import matplotlib as mpl
mpl.use('Agg') #silent mode
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

def plot_work():
    with open('FI-SUM.dat', 'r') as reader:
        lines = reader.readlines()
    lines = [line.strip().split() for line in lines]
    lines = [line for line in lines if line and not line[0] == '#']

    data = np.array(lines, dtype=float)
    steps, works = data[:,0], data[:,4]

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12,8))
    ax.set_title('Energy Change from Free Molecular Dynamics', \
            fontsize=20, fontweight='bold')
    ax.set_xlabel('Step', fontsize=16)
    ax.set_ylabel('Free Energy (Work) / eV', fontsize=16)

    ax.xaxis.set_minor_locator(MultipleLocator(50))
    ax.yaxis.set_minor_locator(MultipleLocator(0.5))

    ax.plot(steps, works, c='b')
    ax.plot(steps, [0]*len(steps), c='k', ls='--') # baseline

    plt.savefig('FI.png')

if __name__ == '__main__':
    plot_work()
