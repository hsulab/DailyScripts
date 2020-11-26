#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse

import numpy as np
norm = np.linalg.norm
dot = np.dot
cross = np.cross

from scipy.interpolate import make_interp_spline, BSpline

import matplotlib as mpl
mpl.use('Agg') #silent mode
from matplotlib import pyplot as plt

from ase.io import read, write

def get_density(density_dict, elements=None):
    """"""
    selected_density = np.zeros(nbins)
    if elements:
        for key, value in density_dict.items():
            if key in elements:
                selected_density += value
    else:
        for key, value in density_dict.items():
            selected_density += value

    return selected_density

def smooth_density(bins, density):
    """"""
    spl = make_interp_spline(bins, density, k=3)
    bins = np.linspace(bins.min(), bins.max(), 300)
    density = spl(bins)

    for i, d in enumerate(density):
        if d < 1e-6:
            density[i] = 0.0

    return bins, density

if __name__ == '__main__':
    atoms = read('CONTCAR-vis', format='vasp')

    # bin setting
    bin_width = 1.0 # bin_center 0.5

    cell = atoms.get_cell()
    nbins = int(np.ceil(cell[2][2]/1.0))
    xyplane_size = norm(cross(cell[0], cell[1]))
    rho_water = 0.997074/18.0152*6.02*1e23*1e-24*18.0152*xyplane_size

    bins = []
    for i in range(nbins):
        bins.append(bin_width/2+bin_width*i)
    bins = np.array(bins)

    # ...
    elements = list(set(atoms.get_chemical_symbols()))
    density_dict = {}
    for e in elements:
        density_dict[e] = np.zeros(nbins)

    atom_masses = atoms.get_masses()
    chemical_symbols = atoms.get_chemical_symbols()
    positions = atoms.get_positions()

    for s, m, p in zip(chemical_symbols, atom_masses, positions):
        for i, bin in enumerate(bins):
            if bin-bin_width/2.0 < p[2] < bin+bin_width/2.0:
                density_dict[s][i] += m #/ xyplane_size

    # plot
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8,8))

    ax.set_title('Mass Density Profile', fontsize=20)
    ax.set_xlabel('Bin Centre / Å')
    ax.set_ylabel('Mass / a.u.')

    # a density profile for water-cluster-surface system
    water_density = get_density(density_dict, ['O','H'])
    cluster_density = get_density(density_dict, ['Pt'])
    surface_density = get_density(density_dict, ['Cd','S'])

    density_dict = {'water': water_density, \
            'cluster': cluster_density,
            'surface': surface_density}

    # smooth
    for name, density in density_dict.items():
        nbins, density = smooth_density(bins, density)
        ax.plot(nbins, density, label=name)

    # baseline
    ax.plot(bins, [0]*len(bins), color='k')

    # waterline
    wbins = [b for b, d in zip(bins, water_density) if d > 1e-6]
    ax.plot(wbins, [rho_water]*len(wbins), 'r_')
    ax.text(wbins[int(len(wbins)*2/3)], rho_water*1.5, r'$\bf{\rho_{H_2O}=1000kg/m^3}$', color='r')

    plt.legend()

    plt.savefig('dp.png')

