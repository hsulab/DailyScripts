#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import numpy as np

import matplotlib as mpl
mpl.use('Agg') #silent mode
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker

import argparse


def plot_reax_training(filename):
    with open(filename, 'r') as reader:
        lines = reader.readlines()
    lines = [line.strip() for line in lines]

    reax_coor, qm_energy, reax_energy = [], [], []
    for line in lines[5:27]:
        line = line.split()
        coor = float(line[2].split('/')[0].split('_')[-1])/100.0
        qm, reax = float(line[4]), float(line[5])
        reax_coor.append(coor)
        qm_energy.append(qm)
        reax_energy.append(reax)

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8,8))
    ax.set_title(r'$C{\equiv}O\ Bond$', fontsize=24)

    ax.set_xlabel('CO Bond Distance / Å', fontsize=20)
    ax.set_ylabel('Relative Energy / kcal/mol', fontsize=20)

    ax.plot(reax_coor, qm_energy, color='r', label='VASP-DFT')
    ax.scatter(reax_coor, qm_energy, marker='s', color='r')

    ax.plot(reax_coor, reax_energy, color='b', label='Lammps-Reax')
    ax.scatter(reax_coor, reax_energy, marker='s', color='b')

    plt.legend()

    plt.savefig('reax_train.png')


if __name__ == '__main__':
    # set argument parser
    parser = argparse.ArgumentParser()

    # add arguments
    parser.add_argument("-f", "--filename", required=True, help="file to plot")

    args = parser.parse_args()

    plot_reax_training(args.filename)
