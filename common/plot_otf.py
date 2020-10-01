#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse

import numpy as np

import matplotlib
matplotlib.use('Agg') #silent mode
import matplotlib.pyplot as plt

MAXLINE = 1000

def iread(fopen):
    results = {}

    results['basic'] = np.nan
    results['corrected'] = np.nan
    results['reference'] = np.nan

    results['called_dft'] = False

    for idx in range(MAXLINE):
        line = fopen.readline()
        if line.startswith('energies'):
            data = np.array((fopen.readline()).strip().split()[1:], dtype=float)
            break
    else:
        data = [np.nan, np.nan, np.nan]
        #raise ValueError('error in reading md step infomation')

    for idx, name in enumerate(['basic','corrected','reference']):
        results[name] = data[idx]

    return results

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-mt', '--metadata', \
            default='hist.dat', help='histogram data file')

    args = parser.parse_args()

    fopen = open('log.txt', 'r')

    energyList = []
    while True:
        line = fopen.readline()
        if line:
            if line.startswith('===== MD Step'):
                results = iread(fopen)
                energyList.append(results)
        else:
            break

    fopen.close()

    nameList = ['basic', 'corrected', 'reference']

    totResults = {}
    totResults['basic'] = {}
    totResults['corrected'] = {}
    totResults['reference'] = {}
    for idx, results in enumerate(energyList):
        for name in nameList:
            if results[name] != np.nan:
                totResults[name][idx] = results[name]

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12,8))
    ax.set_title(
        'On-the-fly MetaDynamics', \
        fontsize=20, fontweight='bold'
    )
    
    ax.set_xlabel('Time / Step', fontsize=16)
    ax.set_ylabel('(Electronic) Free Energy', fontsize=16)

    for name in nameList:
        ax.scatter(
            totResults[name].keys(), 
            totResults[name].values(), 
            label=name
        )
    
    plt.legend()

    plt.savefig('otf.png')
