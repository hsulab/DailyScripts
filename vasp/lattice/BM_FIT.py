#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys

import numpy as np
#import scipy as sp
from scipy.optimize import curve_fit
from scipy.optimize import OptimizeWarning

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


"""
Author: jyxu
Description:
    Brich-Murnaghan Equation (BM Equation)
    E0 - Equilibrium Energy V0 - Equilibrium Volume
    B0 - Bulk Modulus B0' - Derivate of Bulk Modulus
    E(V) = E0 + (9V0B0)/16 * {\
            [(V0/V)^(2/3)-1]^3*B0' + \
            [(V0/V)^(2/3)-1]^2*[6-4(V0/V)^(2/3)]
            }
Notes:
    Currently, only for bulk optimization.
"""


def BM_FIT(data='fitting_data', picname='curve.png'):
    # read system name
    try:
        with open('INCAR', 'r') as reader:
            lines = reader.readlines()
        lines = [line.strip('\n').split() for line in lines]
        sysname = [line[2] for line in lines if len(line) > 0 and line[0] == 'SYSTEM'][0]
        print(sysname)
    except IOError:
        sysname = ''

    # read data
    with open(data, 'r') as reader:
        lines = reader.readlines()
    lines = [line.strip('\n').split() for line in lines]

    dimension = int(lines[0][0][0]) # 2D / 3D
    print('BM Fit for %dD Lattice Optimization.' %dimension)

    x, y = [], []
    for xy in lines[1:]:
        if xy != []:
            x.append(float(xy[0]))
            y.append(float(xy[1]))

    assert len(x) == len(y), 'Inconsitent number of X and Y.'

    # read poscar and calculate volume
    with open('POSCAR', 'r') as reader:
        lines = reader.readlines()
    lines = [line.strip('\n').split() for line in lines]

    a, b, c = [np.array(line, 'float') for line in lines[2:5]]
    volume = np.dot(a, np.cross(b,c)) # experimental volume

    for i, s in enumerate(x):
        x[i] = s**dimension*volume

    # curve fitting
    def BM_FUNCTION(x, E0, V0, B0, B1):
        return E0 + (9.0*V0*B0)/16.0 * (\
                ((V0/x)**(2/3.0)-1)**3.0*B1 + \
                ((V0/x)**(2/3.0)-1)**2.0*(6-4*(V0/x)**(2/3.0))\
                )

    # initial guess is critical for curve fitting
    try:
        #coefs, cov = curve_fit(BM_FUNCTION, x, y)
        initial_guess = [y[int(len(y)/2)], volume, 1., 1.] # E0, V0, B0, B1
        coefs, cov = curve_fit(BM_FUNCTION, x, y, initial_guess)
    except RuntimeWarning or OptimizeWarning:
        print('Initial guess is critical for curve fitting.')
    else:
        print('Initial guess is critical for curve fitting.')

    E0, V0, B0, B1 = coefs
    best_scaling = (V0 / volume)**(1/dimension)

    # save curve-fitting figure
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8,8))
    ax.set_title('%s Lattice Optimization Using BM-Equation' %sysname, \
            fontsize=16, fontweight='bold')

    ax.set_xlabel(r'$\bf{Volume\ /\ Å^3}$', \
            fontsize=16, fontweight='bold')
    ax.set_ylabel('Energy / eV', \
            fontsize=16, fontweight='bold')

    # plot scatter points
    points = ax.scatter(x, y, marker='x', color='r', zorder=20, \
    label='Best Scaling %.2f' %best_scaling)

    # plot fitting curve
    x = np.arange(np.min(x), np.max(x), 0.1)
    y = [BM_FUNCTION(i, E0, V0, B0, B1) for i in x]
    curve = ax.plot(x, y, zorder=10, \
            label='E0 = %.2f V0 = %.2f B0 = %.2f B1 = %.2f ' \
            %(coefs[0], coefs[1], coefs[2], coefs[3]))

    # plot (V0, E0)
    point = ax.scatter(V0, E0, marker='*', color='y', s=50, zorder=30, \
            label='Optimized Lattice')
    
    ax.legend()

    plt.savefig(picname)

if __name__ == '__main__':
    argvs = sys.argv
    nargvs = len(argvs)
    if nargvs == 1:
        BM_FIT()
    elif nargvs == 2:
        data = argvs[1]
        BM_FIT(data)
    elif nargvs == 3:
        data = argvs[1]
        picname = argvs[2]
        BM_FIT(data, picname)
