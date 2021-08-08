#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
plot learning curve 
"""

import argparse

import numpy as np

import matplotlib
matplotlib.use('Agg') #silent mode
import matplotlib.pyplot as plt

from ase.io import read, write 


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--lcurve', 
        default='lcurve.out', 
        help='text file with learning rate and loss'
    )
    parser.add_argument(
        '-s', '--start', type=int,
        default=0, help='text file with learning rate and loss'
    )

    args = parser.parse_args()

    #with open(args.datafile, 'r') as reader:
    #    lines = reader.readlines()

    #lines = [line.strip().split() for line in lines if not line.startswith('#')]

    names = (
        'epoch', 
        'e_trn', 'f_trn', 'e_tst', 'f_tst'
    )
    data = np.loadtxt(args.lcurve)

    fig, axarr = plt.subplots(nrows=2, ncols=1, figsize=(16,12))
    axarr = axarr.flatten()

    LAYOUT =dict(
        left=0.05, right=0.95, 
        bottom=0.05, top=0.95, 
        wspace=0.10,hspace=0.20
    )
    plt.subplots_adjust(**LAYOUT)

    ax = axarr[0]
    ax.set_title(
        'Learning Curve (Test)', 
        fontsize=24, 
        fontweight='bold'
    )

    si = args.start
    
    ax.set_xlabel('Epoch ', fontsize=16)
    ax.set_ylabel('Loss ', fontsize=16)

    ax.plot(data[si:,0], data[si:,3], label='Energy Loss')
    ax.plot(data[si:,0], data[si:,4], label='Force Loss')

    ax.legend(fontsize=24)
    #ax.set_yscale('log')

    ax = axarr[1]
    ax.set_title(
        'Learning Curve (Train)', 
        fontsize=24, 
        fontweight='bold'
    )

    si = args.start
    
    ax.set_xlabel('Batches ', fontsize=16)
    ax.set_ylabel('Loss ', fontsize=16)

    ax.plot(data[si:,0], data[si:,1], label='Energy Loss')
    ax.plot(data[si:,0], data[si:,2], label='Force Loss')

    ax.legend(fontsize=24)
    #ax.set_yscale('log')

    plt.savefig('curve.png')
