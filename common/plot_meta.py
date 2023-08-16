#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.style.use("presentation")

"""Command Line
        plumed sum_hills --stride 100 --hills ../hills.dat --mintozero
"""


def read_fesdat(fesDat):
    """"""
    fopen = open(fesDat, "r")

    # - parse fields...
    fields = fopen.readline().strip().split()
    ncolvar = int((len(fields) - 2) / 2.)

    # - read statistics
    cvDict = {}
    for idx in range(ncolvar):
        cvDict[idx] = {}
        dataLines = []
        for jdx in range(4):
            dataLines.append(fopen.readline())

        cvDict[idx]["min"] = float(dataLines[0].strip().split()[-1])
        cvDict[idx]["max"] = float(dataLines[1].strip().split()[-1])
        cvDict[idx]["nbins"] = int(dataLines[2].strip().split()[-1])
        cvDict[idx]["periodic"] = dataLines[3].strip().split()[-1]

    # - mark start point of the data
    dataLinePoint = fopen.tell()

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
        raise RuntimeError('cannot support number of CV larger than 2')

    freeEnergies = np.array(freeEnergies)

    fopen.close()

    return cvDict, freeEnergies

if __name__ == "__main__":
    """"""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-d", "--datafiles", nargs="*", default=["fes.dat"], 
        help="number of fes.dat files"
    )
    parser.add_argument(
        "--names", nargs="*", default=None, 
        help="names"
    )
    parser.add_argument(
        "-t", "--title", default="MetaDynamics Convergence Analysis", 
        help="title"
    )

    args = parser.parse_args()

    
    datList = args.datafiles
    ndatafiles = len(args.datafiles)
    if args.names is None:
        args.names = [str(i) for i in range(ndatafiles)]
    names = args.names

    kJMol2eV = 1/96.485
    nm2ang = 10.

    # check the CV dimensions
    cvDict, freeEnergies = read_fesdat(datList[0])

    # plot figure(s)
    if len(cvDict.keys()) == 1:
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12,8))
        ax.set_title(args.title)

        ax.set_xlabel("Collective Variable") 
        ax.set_ylabel("Free Energy [eV]") 
        for name, datFile in zip(names, datList):
            cvDict, freeEnergies = read_fesdat(datFile)
            ax.plot(cvDict[0]['points']*nm2ang, freeEnergies*kJMol2eV, label=name)
        ax.legend()
    elif len(cvDict.keys()) == 2:
        for datFile in datList:
            cvDict, free_energies = read_fesdat(datFile)
            free_energies *= kJMol2eV
            cv1, cv2 = cvDict[0]["points"], cvDict[0]["points"]
            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12,8))
            ax.set_xlabel("CV 1")
            ax.set_ylabel("CV 2")
            cn = ax.contour(cv1, cv2, free_energies)
            cntr = ax.contourf(
                cv1, cv2, free_energies, levels=5, #cmap='RdBu'
            )
            fig.colorbar(cntr, ax=ax, label="Free Energy [eV]")
            plt.clabel(cn, inline=True, fontsize=12)
            plt.savefig(os.path.splitext(os.path.basename(datFile))[0]+'.png')
    else:
        pass

    plt.savefig('fes.png')

if __name__ == "__main__":
    ...
