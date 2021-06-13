#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
plot neb curve
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
    parser.add_argument(
        '-df', '--datafile', 
        default='./spline.dat', help='time series files'
    )

    args = parser.parse_args()

    # read data
    data = read_arrays(args.datafile)

    # plot figure
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12,8))
    ax.set_title(
        'NEB Curve',
        fontsize=20, fontweight='bold'
    )
    
    ax.set_xlabel('Reaction Coordinate [Å]', fontsize=16)
    ax.set_ylabel('Potential Energy [eV]', fontsize=16)

    ax.plot(data[:,1], data[:,2])

    plt.savefig('neb.png')
