#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse

import numpy as np

import matplotlib as mpl
mpl.use('Agg') #silent mode
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

"""
Author: Jiayan Xu
Description:
    Plot VASP-MD temperature and energy.
"""

def plot_MD(data_file='MD.dat', \
        fig_title='VASP', fig_format='png'):
    # default figure setting
    timestep = 1.0 # time step 1 fs
    major_interval = 2000
    minor_interval = 1000

    eq_temp = 300 # equilibrium temperature 300 K
    fig_title = r"$\bf{Pt_{13}/CdS(100)}\ $"+"Molecular Dynamics"

    # read data file, the first line is the comment
    with open(data_file, 'r') as reader:
        lines = reader.readlines()
    lines = [line.strip().strip('#').split() for line in lines]

    comments = lines[0]
    data = np.array(lines[1:], dtype=float)

    timesteps = data[:,0]*timestep
    temperatures = data[:,1]
    energies = data[:,2]

    # figure setting
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12,8))
    ax.set_title(fig_title, fontsize=24, fontweight='bold')

    ax.set_xlabel('Time / fs', fontsize=20)
    ax.set_ylabel('Temperature / K', fontsize=20)

    ax.yaxis.set_major_locator(MultipleLocator(20.0))
    ax.yaxis.set_minor_locator(MultipleLocator(5.0))

    t_min, t_max = np.min(temperatures), np.max(temperatures)
    t_min = int(t_min/major_interval) * major_interval
    t_max = (int(t_max/major_interval)+1) * major_interval
    ax.set_xlim(t_min, t_max)

    ax.xaxis.set_major_locator(MultipleLocator(major_interval))
    ax.xaxis.set_minor_locator(MultipleLocator(minor_interval))

    curve_temp, = ax.plot(timesteps, temperatures, color='r') # temperatures
    curve_etemp, = ax.plot(timesteps, [eq_temp]*len(timesteps), \
            color='b', linewidth=2.0) # equilibrium temperature

    ax2 = ax.twinx()
    ax2.set_ylabel('Potential Energy / eV', \
            fontsize=20)

    curve_en, = ax2.plot(timesteps, energies, color='g') # energies

    ax2.yaxis.set_major_locator(MultipleLocator(1.0))
    ax2.yaxis.set_minor_locator(MultipleLocator(0.5))

    plt.legend(handles=[curve_temp, curve_etemp, curve_en], \
            labels=['Temperature', 'Equilibrium Temperature', 'Potential Energy'])

    plt.savefig('MD.'+fig_format)

if __name__ == '__main__':
    #parser = argparse.ArgumentParser()

    #parser.add_argument('-f', '--datafile', \
    #        default='PDOS_USER.dat', help='PDOS Data File')
    #parser.add_argument('-o', '--out', \
    #        default='DOS.png', help='DOS Figure')

    #args = parser.parse_args()

    plot_MD()
