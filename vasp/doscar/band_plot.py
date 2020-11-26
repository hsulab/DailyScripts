#!/usr/bin/env python3
# -*- coding:utf-8 -*-

import os
import sys

import numpy as np
import matplotlib as mpl
mpl.use('Agg') #silent mode
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker


if __name__ == '__main__':
    # read k-point labels
    group_labels, x = [], []
    
    with open("KLABELS", 'r') as reader:
        lines = reader.readlines()[1:]

    for i in lines:
        s = i.encode('utf-8')#.decode('latin-1')
        if len(s.split()) == 2 and not s.decode('utf-8','ignore').startswith("*"):
            group_labels.append(s.decode('utf-8','ignore').split()[0])
            x.append(float(s.split()[1]))

    # read band gap
    with open("BAND_GAP", 'r') as reader:
        lines = reader.readlines()
    lines = [line.strip('\n').split() for line in lines]
    band_gap = float(lines[2][-1])
    
    #hsp=np.loadtxt(open("KLABELS", encoding='utf8'), dtype=np.string_,skiprows=1,usecols = (0,1))
    for index in range(len(group_labels)):
        if group_labels[index] == "GAMMA":
            group_labels[index] = u"Γ"

    # read INCAR SYSTEM
    with open('INCAR', 'r') as reader:
        lines = reader.readlines()
    lines = [line.strip().split('=') for line in lines]
    system_name = [line[1].split('#')[0] for line in lines if line[0].strip(' ') == 'SYSTEM'][0]
    
    # plot bands
    fig, axe = plt.subplots(nrows=1, ncols=1, figsize=(8,8))
    axe.set_title('Band Structure of %s \nBand Gap is %.3f eV' \
            %(system_name, band_gap), fontsize=20, fontweight='bold')

    '''
    df_up = np.loadtxt("REFORMATTED_BAND_UP.dat", dtype=np.float64)
    df_dw = np.loadtxt("REFORMATTED_BAND_DW.dat", dtype=np.float64)

    df = df_up
    axe.plot(df[:,0], df[:,1:], linewidth=1.0, color='red')
    
    df = df_dw
    axe.plot(df[:,0], df[:,1:], linewidth=1.0, color='blue')
    '''
    if os.path.exists('REFORMATTED_BAND_UP.dat') and \
            os.path.exists('REFORMATTED_BAND_DW.dat'):
        datfiles = ['REFORMATTED_BAND_UP.dat', 'REFORMATTED_BAND_DW.dat']
        colors = ['red', 'blue']
    elif os.path.exists('REFORMATTED_BAND.dat'):
        datfiles = ['REFORMATTED_BAND.dat']
        colors = ['blue']
    else:
        raise ValueError('No such data files.')

    for dat, color in zip(datfiles, colors):
        df = np.loadtxt(dat, dtype=np.float64)
        axe.plot(df[:,0], df[:,1:], linewidth=2.0, color=color)

    # plot labels and 
    font = {'color'  : 'black',  
            'weight' : 'bold',  
            'size'   : 14,  
            }  
    
    axe.set_xlabel(r'$\it{\bf{\mathbf{k}-points}}$', fontdict=font)
    axe.set_ylabel(r'$\it{\bf{{E}-E_{f}\ (eV)}}$', fontdict=font)

    axe.set_xticks(x)

    axe.set_xticklabels(group_labels, rotation=0, fontdict=font)

    # plot high symmetry k-point
    axe.axhline(y=0, xmin=0, xmax=1,linestyle= '--',linewidth=1.5,color='0.5')
    for i in x[1:-1]:
    	axe.axvline(x=i, ymin=0, ymax=1,linestyle= '--',linewidth=1.5,color='0.5')
    axe.set_xlim((x[0], x[-1]))
    
    #plt.show()
    plt.savefig('band.png')
