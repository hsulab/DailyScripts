#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" Plot metadata
"""

import os
import argparse

import numpy as np

import matplotlib
matplotlib.use('Agg') #silent mode
import matplotlib.pyplot as plt

def read_metadata(infile):
    """read metadata"""
    with open(infile, 'r') as reader:
        lines = reader.readlines()
    lines = [
        line.strip().split() 
        for line in lines if not line.strip().startswith('#')
    ]
    datfiles = [line[0] for line in lines]
    coords = [float(line[1]) for line in lines]

    return datfiles, coords

def read_usdat(datfile, nEquil):
    """ read usdat
    """
    with open(datfile, 'r') as reader:
        lines = reader.readlines()
    lines = [
        line.strip().split() for line in lines 
        if not line.startswith('#')
    ]

    data = np.array(lines, dtype='float')[nEquil:, 1]

    nsamples = len(data)
    smin, smax = np.min(data), np.max(data)
    if np.fabs(smin-smax) > 0.02:
        winfo = 'WARN!!!'
    else:
        winfo = ""
    info = ("{:<20s}{:>8d}{:>8.4f}{:>8.4f}{:>10s}").format(datfile, nsamples, smin, smax, winfo)
    print(info)

    return data

def plot_bins(metadata, nBins, nEquil):
    datfiles, coords = read_metadata(metadata)

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12,8))
    ax.set_title('Umbrella Sampling Overlap Analysis', \
            fontsize=20, fontweight='bold')

    ax.set_xlabel('Collective Variable', fontsize=16)
    ax.set_ylabel('Number', fontsize=16)

    frames = []
    for datfile in datfiles:
        if os.path.exists(datfile):
            data = read_usdat(datfile, nEquil)
            n, bins, patches = ax.hist(data, bins=nBins, alpha=0.45)
            bins = np.array(bins)
            centres = (bins[1:] + bins[:-1])/2.
            ax.plot(centres,n,'k-',lw=2)
            frames.append(data)
        else:
            print('Not exists %s' %datfile)

    frames = np.array(frames)

    content = '# hist\n'
    for line in frames.T:
        content += ('{:<5.3f}  '*len(line)+'\n').format(*list(line))
    with open('hist.dat','w') as writer:
        writer.write(content)

    ax.scatter(coords, len(coords)*[0], marker='*', color='k')

    plt.savefig('us-bins.png')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-mt', '--metadata', \
            default='hist.dat', help='histogram data file')
    parser.add_argument('-nb', '--numbins', type=int,\
            default=10, help='number of bins')
    parser.add_argument('-ne', '--nequil', type=int,\
            default=300, help='number of equilibrium steps')

    args = parser.parse_args()

    plot_bins(args.metadata, args.numbins, args.nequil)
