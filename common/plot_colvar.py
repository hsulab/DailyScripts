#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
plot time versus colvar
"""

import argparse

import numpy as np

import matplotlib
matplotlib.use('Agg') #silent mode
import matplotlib.pyplot as plt

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-df', '--datafile', \
            default='bias', help='time series files')

    args = parser.parse_args()

    with open(args.datafile, 'r') as reader:
        lines = reader.readlines()

    lines = [line.strip().split() for line in lines if not line.startswith('#')]

    data = np.array(lines, dtype=float)

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12,8))
    ax.set_title(
        'Collective Variable', \
        fontsize=20, fontweight='bold'
    )
    
    ax.set_xlabel('Time / Step', fontsize=16)
    ax.set_ylabel('Collective Variable', fontsize=16)

    plt.plot(data[:,0], data[:,1])

    plt.savefig('colvar.png')
