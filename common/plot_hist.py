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

TITLE = r"Histogram"
XLABEL = r"Collective Variable"
YLABEL = "Number"


def plot_hist(datfile='hist.dat',pic='hist.png',cols=[]):
    # read hist data
    with open(datfile, 'r') as reader:
        lines = reader.readlines()
    lines = [line.strip().split() for line in lines \
            if not line.strip().startswith('#')]
    data = np.array(lines, dtype='float')

    data = data.T

    # get columns
    if cols:
        frames = []
        for i in cols:
            frames.append(data[i])
    else:
        frames = data

    # figure
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=FIGSIZE)

    # labels
    ax.set_title(r"Histogram",TFONT)
    ax.set_xlabel(XLABEL, LFONT)
    ax.set_ylabel("Number", LFONT)

    # plot histogram
    for frame in frames:
        ax.hist(frame,bins=10,alpha=0.45)

    fig.tight_layout()

    plt.savefig(pic)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-d', '--dat', \
            default='hist.dat', help='histogram data file')
    parser.add_argument('-p', '--pic', \
            default='hist.png', help='output picture file')
    parser.add_argument('-c', '--col', \
            nargs='*', type=int, help='index of column in Dat (start from ZERO)')

    args = parser.parse_args()

    plot_hist(args.dat,args.pic,args.col)
