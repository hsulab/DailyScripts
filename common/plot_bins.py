#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" Plot metadata
"""

import argparse

import numpy as np

import matplotlib
matplotlib.use('Agg') #silent mode
import matplotlib.pyplot as plt

def read_metadata(infile):
    """read metadata"""
    with open(infile, 'r') as reader:
        lines = reader.readlines()
    lines = [line.strip().split() \
            for line in lines if not line.strip().startswith('#')]
    datfiles = [line[0] for line in lines]
    coords = [float(line[1]) for line in lines]


    return datfiles, coords

def read_usdat(datfile):
    """read usdat
    """
    with open(datfile, 'r') as reader:
        lines = reader.readlines()
    lines = [line.strip().split() for line in lines \
            if not line.startswith('#')]

    data = np.array(lines, dtype='float')[:, 1]

    nsamples = len(data)
    smin, smax = np.min(data), np.max(data)
    if np.fabs(smin-smax) > 0.02:
        winfo = 'WARN!!!'
    else:
        winfo = ""
    info = ("{:<20s}{:>8d}{:>8.4f}{:>8.4f}{:>10s}").format(datfile, nsamples, smin, smax, winfo)
    print(info)

    return data

def plot_bins(metadata):
    datfiles, coords = read_metadata(metadata)

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
    parser = argparse.ArgumentParser()
    parser.add_argument('-mt', '--metadata', \
            default='hist.dat', help='histogram data file')

    args = parser.parse_args()

    plot_bins(args.metadata)
