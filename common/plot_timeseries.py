#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
plot time versus colvar
the first column (0) must be the time if required
"""

import argparse

import numpy as np

import matplotlib
matplotlib.use('Agg') #silent mode
import matplotlib.pyplot as plt

def read_arrays(datafile):
    with open(datafile, 'r') as reader:
        lines = reader.readlines()

    lines = [line.strip().split() for line in lines if not line.startswith('#')]

    data = np.array(lines, dtype=float)

    return data

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-df', '--datafile', \
            default='bias', help='time series files')
    parser.add_argument('-c', '--column', nargs='*', type=int, \
            default=[1], help='time series files')

    args = parser.parse_args()

    # read data
    data = read_arrays(args.datafile)

    # plot figure
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12,8))
    ax.set_title(
        'Collective Variable Time Series', \
        fontsize=20, fontweight='bold'
    )
    
    ax.set_xlabel('Time / Step', fontsize=16)
    ax.set_ylabel('Collective Variable', fontsize=16)

    for col in args.column:
        ax.plot(range(len(data[:,col])), data[:,col])

    plt.savefig('series.png')
