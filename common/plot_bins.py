#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse

import numpy as np

import matplotlib as mpl
mpl.use('Agg') #silent mode
from matplotlib import pyplot as plt

def read_metadata(infile):
    """"""
    with open(infile, 'r') as reader:
        lines = reader.readlines()
    lines = [line.strip().split() \
            for line in lines if not line.strip().startswith('#')]
    datfiles = [line[0] for line in lines]
    coords = [float(line[1]) for line in lines]

    return datfiles, coords

def read_usdat(datfile):
    with open(datfile, 'r') as reader:
        lines = reader.readlines()
    lines = [line.strip().split() for line in lines \
            if not line.startswith('#')]

    data = np.array(lines, dtype='float')[:,1]

    return data

def plot_bins():
    datfiles, coords = read_metadata('METADATA')

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12,8))
    ax.set_title('Umbrella Sampling Overlap Analysis', \
            fontsize=20, fontweight='bold')

    ax.set_xlabel('Collective Variable', fontsize=16)
    ax.set_ylabel('Number', fontsize=16)

    frames = []
    for datfile in datfiles:
        data = read_usdat(datfile)
        n, bins, patches = ax.hist(data, bins=10, alpha=0.45)
        bins = np.array(bins)
        centres = (bins[1:] + bins[:-1])/2.
        #print(n)
        #print(bins)
        ax.plot(centres,n,'k-',lw=2)
        #print(patches[0].get_color())

        frames.append(data)
    frames = np.array(frames)

    content = '# hist\n'
    for line in frames.T:
        content += ('{:<5.3f}  '*len(line)+'\n').format(*list(line))
    with open('hist.dat','w') as writer:
        writer.write(content)

    ax.scatter(coords, len(coords)*[0], marker='*', color='k')

    plt.savefig('us-bins.png')

if __name__ == '__main__':
    plot_bins()
