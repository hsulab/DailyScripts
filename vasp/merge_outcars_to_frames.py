#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse

from ase.io import read, write
from ase.io.vasp import read_vasp_out

if __name__ == '__main__':
    # args 
    parser = argparse.ArgumentParser()
    parser.add_argument('-mt', '--metadata', \
            default='hist.dat', help='histogram data file')
    parser.add_argument('-nb', '--numbins', type=int,\
            default=10, help='number of bins')
    parser.add_argument('-ne', '--nequil', type=int,\
            default=300, help='number of equilibrium steps')

    args = parser.parse_args()

    # data 
    data_files = ['mt-x1/OUTCAR', 'mt-x11/OUTCAR', 'mt-x12/OUTCAR']

    frames = []
    for file_name in data_files:
        frames.extend(read_vasp_out(data_files[0], ':'))

    # set steps
    for idx, atoms in enumerate(frames):
        atoms.info['step'] = idx

    write('total_outcars.xyz', frames)
