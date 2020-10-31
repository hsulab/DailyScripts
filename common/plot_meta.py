#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse

import numpy as np

import matplotlib
matplotlib.use('Agg') #silent mode
import matplotlib.pyplot as plt


def read_fesdat(fesDat):
    #with open(fesDat, 'r') as reader:
    #    lines = reader.readlines()

    #data = np.array(
    #    [
    #        line.strip().split() 
    #        for line in lines if not line.startswith('#!')
    #    ], dtype=float
    #)

    fopen = open(fesDat, 'r')

    fields = fopen.readline().strip().split()
    ncolvar = int((len(fields) - 2) / 2.)

    cvDict = {}
    for idx in range(ncolvar):
        cvDict[idx] = {}
        dataLines = []
        for jdx in range(4):
            dataLines.append(fopen.readline())
        
        cvDict[idx]['min'] = float(dataLines[0].strip().split()[-1])
        cvDict[idx]['max'] = float(dataLines[1].strip().split()[-1])
        cvDict[idx]['nbins'] = int(dataLines[2].strip().split()[-1])
        cvDict[idx]['periodic'] = dataLines[3].strip().split()[-1]

        print(cvDict[idx])

    dataLinePoint = fopen.tell()

    #print(fopen.readline())

    # get grids
    lastNumBins = 0
    for idx in range(ncolvar):
        fopen.seek(dataLinePoint)
        dataLists = []
        for jdx in range(cvDict[idx]['nbins']):
            dataLists.append(fopen.readline().strip().split())
            for kdx in range(lastNumBins):
                line = fopen.readline()
        dataList = [data[idx] for data in dataLists]
        cvDict[idx]['points'] = np.array(dataList, dtype=float)
        if lastNumBins == 0:
            lastNumBins = cvDict[idx]['nbins']
        else:
            lastNumBins = lastNumBins*cvDict[idx]['nbins']

    # free energy
    # TODOL: only 2 CV for now, maybe extended later
    freeEnergies = []
    if ncolvar == 1:
        fopen.seek(dataLinePoint)
        curList = []
        for idx in range(cvDict[0]['nbins']):
            curList.append(float(fopen.readline().strip().split()[1]))
        freeEnergies.extend(curList)
    elif ncolvar == 2:
        fopen.seek(dataLinePoint)
        for idx in range(cvDict[1]['nbins']):
            curList = []
            for jdx in range(cvDict[0]['nbins']):
                curList.append(float(fopen.readline().strip().split()[2]))
            freeEnergies.append(curList)
            emptyLine = fopen.readline()
    else:
        raise ValueError('cannot support number of CV larger than 2')

    freeEnergies = np.array(freeEnergies)

    fopen.close()

    return cvDict, freeEnergies

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-nf', '--numfiles', type=int,\
            default=1, help='number of fes.dat files')

    args = parser.parse_args()

    datList = []
    if args.numfiles == 1:
        datList.append('fes.dat')
    else:
        for idx in range(args.numfiles):
            datList.append('fes_'+str(idx)+'.dat')

    # check the CV dimensions
    cvDict, freeEnergies = read_fesdat(datList[0])

    # plot figure(s)
    if len(cvDict.keys()) == 1:
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12,8))
        ax.set_title('MetaDynamics Convergence Analysis', \
                fontsize=20, fontweight='bold')

        ax.set_xlabel('Collective Variable', fontsize=16)
        ax.set_ylabel('Free Energy', fontsize=16)

        for datFile in datList:
            cvDict, freeEnergies = read_fesdat(datFile)
            ax.plot(cvDict[0]['points']*10, freeEnergies, label=datFile)
        ax.legend()
    elif len(cvDict.keys()) == 2:
        for datFile in datList:
            cvDict, freeEnergies = read_fesdat(datFile)
            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12,8))
            cn = ax.contour(cvDict[0]['points'], cvDict[1]['points'], freeEnergies) 
                #levels=14, linewidth=0.5, color='k')
            cntr = ax.contourf(cvDict[0]['points'], cvDict[1]['points'], freeEnergies, 
                levels=12, cmap='RdBu')
            fig.colorbar(cntr, ax=ax)
            plt.clabel(cn, inline=True, fontsize=12)
            plt.savefig(os.path.splitext(os.path.basename(datFile))[0]+'.png')
    else:
        pass

    plt.savefig('fes.png')
