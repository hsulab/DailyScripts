#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse

import numpy as np
import pandas as pd

import pyblock as pb

def block_analysis(datfile='BM-1.dat', ncol=2, ndrop=1000):
    with open(datfile, 'r') as reader:
        lines = reader.readlines()
    lines = [line.strip().split() for line in lines if not line.startswith('#')]

    print('Successfully read blue moon data from %s.' %datfile)

    # 
    data = np.array(lines, dtype=float)

    data = pd.Series(data[ndrop:,ncol]) # gradients

    (data_length, reblock_data, covariance) = pb.pd_utils.reblock(data)
    reblock_data.to_csv('RB.csv')

    pb.plot.plot_reblocking(reblock_data, 'RB.png', False)
    print('Successfully carry out block analysis.')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-f', '--file', help='Bluemoon File')
    parser.add_argument('-nd', '--ndrop', type=int, help='Number of Drops')
    parser.add_argument('-nc', '--ncol', type=int, default=3, \
            help='Number of Drops')

    args = parser.parse_args()

    block_analysis(args.file, args.ncol-1, args.ndrop)
