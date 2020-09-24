#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse

import numpy as np

import matplotlib as mpl
mpl.use('Agg') #silent mode
from matplotlib import pyplot as plt

from matplotlib import cm


def read_whamdat(whamdat):
    nbins_x, nbins_y = 84, 8

    fopen = open(whamdat, 'r')
    line = fopen.readline() # first line of file
    data = []
    for i in range(nbins_x):
        cur_data = []
        for j in range(nbins_y):
            line = fopen.readline()
            if line:
                cur_data.append(line.split()[:3])
            else:
                print('End of file.')
                break
        line = fopen.readline() # blank line
        data.append(cur_data)
    fopen.close()

    # data
    data = np.array(data, dtype=float)

    X, Y, Z = [], [], []
    for cur_data in data:
        X.append(cur_data[0][0])
        cur_z = []
        for cur_row in cur_data:
            cur_z.append(cur_row[2])
        Z.append(cur_z)
    Y = data[0][:,1]

    X = np.array(X)
    Z = np.array(Z).T

    mZ = np.ma.masked_equal(Z,np.inf)
    #Z = np.ma.array(Z, mask=mz.mask)

    return X, Y, mZ


def plot_wham():
    # read data
    x, y, e = read_whamdat('us-new.csv')
    e = e / 23.061 # kcal/mol2eV

    # plot figure
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12,8))
    ax.set_title('Path Collective Variable Umbrella Sampling', \
            fontsize=20, fontweight='bold')
    ax.set_xlabel('s (progress along reaction path)', fontsize=16)
    ax.set_ylabel('z (distance away from path)', fontsize=16)

    levels = np.arange(np.min(e), np.max(e), 0.5)

    norm = cm.colors.Normalize(vmax=abs(e).max(), vmin=-abs(e).max())
    cmap = cm.PRGn

    #ax.contourf(x, y, e, levels, norm=norm, \
    #        cmap=cm.get_cmap(cmap, len(levels)-1))
    cs = ax.contourf(x, y, e, 10)

    cl = ax.contour(x, y, e, cs.levels, colors='k')
    for c in cl.collections:
        c.set_linestyle('solid')

    # contour label
    ax.clabel(cl, fmt='%2.1f', colors='k', fontsize=14)

    # color bar
    cbar = fig.colorbar(cs)
    cbar.ax.set_ylabel('Free Energy (kcal / mol)', fontsize=16)

    # find MFEP
    paths = []
    for i, e_row in enumerate(e.T):
        p = 0
        emin = 1e8
        for j, cur_e in enumerate(e_row):
            if cur_e < emin:
                p = j
        paths.append([x[i],y[p]])
    paths = np.array(paths)

    ax.plot(paths[:,0], paths[:,1], linestyle='solid', color='k', linewidth=4)

    plt.savefig('f2d.png')

if __name__ == '__main__':
    '''
    parser = argparse.ArgumentParser()

    parser.add_argument('-f', '--datafile', \
            default='PDOS_USER.dat', help='PDOS Data File')
    parser.add_argument('-o', '--out', \
            default='DOS.png', help='DOS Figure')

    args = parser.parse_args()
    '''

    #read_whamdat('us.csv')
    plot_wham()
