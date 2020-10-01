#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse

import numpy as np

import matplotlib
matplotlib.use('Agg') #silent mode
import matplotlib.pyplot as plt

def read_fesdat(fesDat):

    with open(fesDat, 'r') as reader:
        lines = reader.readlines()

    data = np.array(
        [
            line.strip().split() 
            for line in lines if not line.startswith('#!')
        ], dtype=float
    )

    return data

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-nf', '--numfiles', type=int,\
            default=1, help='number of fes.dat files')

    args = parser.parse_args()

    datList = []
    if args.numfiles == 1:
        datList.append('fes.dat')
    else:
        for idx in range(args.numfiles):
            datList.append('fes_'+str(idx)+'.dat')
    print(datList)

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12,8))
    ax.set_title('MetaDynamics Convergence Analysis', \
            fontsize=20, fontweight='bold')

    ax.set_xlabel('Collective Variable', fontsize=16)
    ax.set_ylabel('Free Energy', fontsize=16)

    for datFile in datList:
        data = read_fesdat(datFile)
        ax.plot(data[:,0]*10, data[:,1], label=datFile)

    ax.legend()

    plt.savefig('fes.png')
