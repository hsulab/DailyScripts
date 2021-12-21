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
plt.style.use("presentation")

from ase.io import read, write

def get_density(density_dict, nbins, elements=None):
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

def calc_density(atoms, merged=False):
    # bin setting
    bin_width = 1.0 # bin_center 0.5

    # calculate xy area
    cell = atoms.get_cell()
    nbins = int(np.ceil(cell[2][2]/bin_width))
    xyplane_size = norm(cross(cell[0], cell[1]))

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
        for i, bin_pos in enumerate(bins):
            if bin_pos-bin_width/2.0 < p[2] < bin_pos+bin_width/2.0: # use z-axis
                density_dict[s][i] += m #/ xyplane_size

    if merged:
        density_data = get_density(density_dict, nbins)
    else:
        density_data = density_dict

    return density_data, bins

def plot_density(density_data):
    # plot
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(16,12))

    plt.suptitle('Mass Density Profile')
    ax.set_xlabel('Bin Centre [Å]')
    ax.set_ylabel('Mass [a.u.]')

    # a density profile for water-cluster-surface system
    #water_density = get_density(density_dict, ['O','H'])
    #cluster_density = get_density(density_dict, ['Pt'])
    #surface_density = get_density(density_dict, ['Cd','S'])

    #density_dict = {
    #    'water': water_density, 
    #    'cluster': cluster_density,
    #    'surface': surface_density
    #}

    #density_dict = {
    #    "all": get_density(density_dict, nbins)
    #}

    # smooth
    #for name, density in density_dict.items():
    #    nbins, density = smooth_density(bins, density)
    #    ax.plot(nbins, density, label=name)

    nbins, density = smooth_density(bins, density_data)
    ax.plot(nbins, density, color="k", label="avg")

    nbins, density = smooth_density(bins, density_max)
    ax.plot(nbins, density, ls="--", color="grey", label="max")

    nbins, density = smooth_density(bins, density_min)
    ax.plot(nbins, density, ls="-.", color="grey", label="min")

    # baseline
    ax.plot(bins, [0]*len(bins), color='k')

    # waterline
    #  solvent density [g/cm3] / molecule weight [g/mol] * [mol] * [A2cm]^3 * [g/mol] * [A^2]
    #rho_solvent = 0.997074*1e-24*6.02*1e23*xyplane_size

    #wbins = [b for b, d in zip(bins, water_density) if d > 1e-6]
    #ax.plot(wbins, [rho_water]*len(wbins), 'r_')
    #ax.text(wbins[int(len(wbins)*2/3)], rho_water*1.5, r'$\bf{\rho_{H_2O}=1000kg/m^3}$', color='r')

    plt.legend()

    plt.savefig('dp.png')


    return

def plot_density_file(density_file="./density.dat"):
    #np.

    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "structure"
    )
    args = parser.parse_args()

    frames = read(args.structure, ":")
    
    all_data = []
    for i, atoms in enumerate(frames):
        density_data, bins = calc_density(atoms, True)
        #print(density_data)
        # print(density_data)
        all_data.append(density_data)
    all_data = np.array(all_data)
    density_max = np.max(all_data, axis=0) # max along time
    density_min = np.min(all_data, axis=0) # min along time
    density_avg = np.mean(all_data, axis=0)

    all_data = np.vstack([all_data, density_min, density_avg, density_max])

    np.savetxt("density.dat", all_data, fmt="%8.4f")

    plot_density(density_avg)
