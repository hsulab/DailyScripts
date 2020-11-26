#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse

import numpy as np

import matplotlib as mpl
mpl.use('Agg') #silent mode
from matplotlib import pyplot as plt

def plot_feg(datfile, fname, pname):
    with open(datfile, 'r') as reader:
        lines = reader.readlines()

    lines_new = []
    for line in lines:
        if line.startswith('#'):
            print(line)
        else:
            line = line.strip().split()
            lines_new.append(line)
    lines = lines_new
    data = np.array(lines, dtype=float)

    steps = data[:,0]
    coords = data[:,1]
    grads = data[:,2]

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12,8))
    ax.set_title(fname+' Free Energy Gradient Change', \
            fontsize=20, fontweight='bold')
    ax.set_xlabel('Molecular Dynamics Step', fontsize=16)
    ax.set_ylabel('Free Energy Gradient / (dA/dξ)', fontsize=16)

    a, = ax.plot(steps[1:], grads[1:], color='b', label='Gradient')

    ax2 = ax.twinx()
    b, = ax2.plot(steps, coords, color='r', label='Coordinate')
    ax2.set_ylabel('Reaction Coordinate / a.u.', fontsize=16)

    plt.legend(handles=[a,b])

    plt.savefig(pname)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-d', '--datfile', nargs='?',\
            default='ffff.dat', help='PDOS Data File')
    parser.add_argument('-f', '--fname', nargs='?',\
            default='ffff', help='PDOS Data File')
    parser.add_argument('-p', '--pname', nargs='?',\
            default='ffff.png', help='DOS Figure')

    args = parser.parse_args()

    plot_feg(args.datfile, args.fname, args.pname)
