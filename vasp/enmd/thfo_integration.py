#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse

import numpy as np

from scipy import integrate

import argparse

import numpy as np

import matplotlib as mpl
mpl.use('Agg') #silent mode
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

def thermo_integration(datfile='THFO-1.dat', nframe=500, intv=0.1, region=None):
    # read BM.dat file
    with open(datfile, 'r') as reader:
        lines = reader.readlines()
    # lines = [line.strip().split() for line in lines if not line.startswith('#')]

    lines_new = []
    for line in lines:
        if line.startswith('#'):
            print(line)
        else:
            line = line.strip().split()
            lines_new.append(line)
    lines = lines_new
    data = np.array(lines, dtype=float)

    print('Successfully read ThermoDynamic data from %s.' %datfile)

    coords = data[:nframe,1] # collective variable
    gradients = data[:nframe,2] # free energy gradients

    if coords[0] - coords[-1] > 0.0:
        findmax = 1.0
    else:
        findmax = -1.0

    # bins
    rc_min, rc_max = np.min(coords), np.max(coords)
    rc_min = np.floor(rc_min*10.0)/10.0
    rc_max = np.ceil(rc_max*10.0)/10.0

    bins = np.arange(rc_min,rc_max,intv) # start of each bin

    nintvs = int(np.ceil((rc_max - rc_min)/intv))
    print(rc_min, rc_max, nintvs)

    his = [[] for i in range(nintvs)] # accumulated gradients

    for i, start in enumerate(bins):
        for coord, grad in zip(coords, gradients):
            if start < coord < start + intv:
                his[i].append(grad)

    grad_avg = np.zeros(nintvs)
    std_errs = np.zeros(nintvs)
    for i, h in enumerate(his):
        if h:
            grad_avg[i]=np.mean(h)
            std_errs[i]=np.std(h)
        else:
            grad_avg[i]=0
            std_errs[i]=0

    coord_avg = np.zeros(nintvs)
    for i in range(nintvs):
        coord_avg[i] = bins[i] + intv/2.0

    # remove bin without samples
    coord_new, grad_new, std_new = [], [], []
    for coord, grad, err in zip(coord_avg, grad_avg, std_errs):
        if grad != 0:
            coord_new.append(coord)
            grad_new.append(grad)
            std_new.append(err)
    coord_avg, grad_avg, std_errs = coord_new, grad_new, std_new
    
    # set integrate and plot region
    if region:
        low, high = region[0], region[1]
        coord_new, grad_new = [], []
        for coord, grad in zip(coord_avg,grad_avg):
            if low < coord < high:
                coord_new.append(coord)
                grad_new.append(grad)
        coord_avg, grad_avg = coord_new, grad_new
        print('Set reactive coordinates between %.2f and %.2f.' %(low, high))
    # print(coord_avg)
    # print(grad_avg)

    feind, femark = 0, 0.0
    free_energies = []
    for i in range(len(grad_avg)):
        coord, grad = coord_avg[:i+1], grad_avg[:i+1]
        v = integrate.trapz(grad, coord)
        free_energies.append(v)
        if findmax*v >= femark:
            femark = findmax*v
            feind = i

    # print(free_energies)
    print('Successfully carry out the thermointegration.')

    # plot the figure
    fig, ax = plt.subplots(1, 1, figsize=(12,8))
    ax.set_title('Thermodynamic Force and Integration', fontsize=20, fontweight='bold')

    ax.set_xlabel('Reactive Coordinates / a.u.', fontsize=16)
    ax.set_ylabel('Thermo Force (dA/dξ)', fontsize=16)
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))

    ymin = np.floor(np.min(grad_avg)*10.0)/10.0
    ymax = np.ceil(np.max(grad_avg)*10.0)/10.0

    ymin2 = np.floor(np.min(free_energies)*10.0)/10.0
    ymax2 = np.ceil(np.max(free_energies)*10.0)/10.0

    ymin = np.min([ymin, ymin2])
    ymax = np.max([ymax, ymax2])

    ax.set_ylim(ymin, ymax)
    ax.yaxis.set_minor_locator(MultipleLocator(0.1))

    # ax.scatter(coord_avg, grad_avg, marker='x', c='r')
    ax.plot(coord_avg, [0]*len(coord_avg), ls='--', c='k')
    a, = ax.plot(coord_avg, grad_avg, ls='-', c='b', label='Thermo Force')

    for c, g, h in zip(coord_avg, grad_avg, his):
        ax.text(c*1.01,g*1.01,'%d' %len(h), fontweight='bold')

    ax.errorbar(coord_avg, grad_avg, std_errs, \
            marker='s', mfc='g', mec='g', ecolor='g')

    ax2 = ax.twinx()
    ax2.set_ylabel('Free Energy / eV', fontsize=16)
    ax2.set_ylim(ymin, ymax)
    ax2.yaxis.set_minor_locator(MultipleLocator(0.1))

    b, = ax2.plot(coord_avg, free_energies, color='r', label='Free Energy / eV')
    ax2.scatter(coord_avg[feind], findmax*femark, marker='*', c='y') # mark IS or TS
    ax2.text(coord_avg[feind]*1.01,findmax*femark*1.01,'%.2f eV' %(findmax*femark), fontweight='bold')

    # ax.set_xticks(np.arange(0, np.max(steps), 0.5))
    # ax2.plot(coords, [0]*len(coords), c='k', ls='--')

    plt.legend(handles=[a, b], loc='best')

    plt.savefig('THFO.png')

    print('Successfully plot the thermointegration process to THFO.png.')

    # write free energies to file
    content = '# step, gradient, energy\n'
    for step, force, energy in zip(coord_avg, grad_avg, free_energies):
        content += '{:<12.4f}{:<12.4f}{:<12.4f}\n'.format(step, force, energy)

    with open('THFOI.dat', 'w') as writer:
        writer.write(content)

    print('Successfully write thermointegration data into THFOI.dat.')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file', required=True, \
            help='BM datfile')
    parser.add_argument('-nf', '--nframes', required=True, \
            type=int, help='Number of Frames')
    parser.add_argument('-r', '--region', nargs='*', default=None, \
            type=float, help='Reactive Coordinates Range')
    parser.add_argument('-bw', '--binwidth', nargs='?', default=0.1, \
            type=float, help='Bin Width')
    args = parser.parse_args()

    thermo_integration(args.file, args.nframes, args.binwidth, args.region)
    # thermo_integration('THFO-1.dat', 500)
