#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse

import numpy as np

import matplotlib as mpl
mpl.use('Agg') #silent mode
from matplotlib import pyplot as plt

"""
Author: Jiayan Xu
"""

def read(filename='TFILOG', ncons=2, tolstep=3):
    # get file object
    reader = open(filename, 'r')
        
    nstep, flag = 0, True
    tot_rcs, tot_thfos = [], []
    cur_rcs, cur_thfos = [], []
    while flag:
        # find MD step
        line = reader.readline().strip()
        if line.startswith('--------------------    MD STEP'):
            # count and reset data
            nstep += 1
            cur_rcs, cur_thfos = [], []
            irbm = False # If Read Blue Moon

        # read section
        if line.startswith('-> REACTIVE'):
            line = reader.readline().strip() # blank line
            for i in range(ncons):
                line = reader.readline().strip()
                val = float(line.split()[1])
                cur_rcs.append(val)
            tot_rcs.append(cur_rcs)

        if line.startswith('-> FREE ENERGY'):
            line = reader.readline().strip() # blank line
            line = reader.readline().strip() # RC, FEG, FEG1, FEG2
            for i in range(ncons):
                line = reader.readline().strip()
                val = float(line.split()[2])
                cur_thfos.append(val)
            tot_thfos.append(cur_thfos)
            irbm = True

        # check step
        if nstep == tolstep and irbm == True:
            flag = False

    reader.close()
    print('Successfully read CONSTRAINTS and LAMBDAS in TFILOG.')

    # check output
    tot_rcs = np.array(tot_rcs)
    tot_thfos = np.array(tot_thfos)

    # print(tot_rcs)

    # write constraints and lambdas to THFO.dat
    for i in range(ncons):
        content = ('{:<8s}'+'{:<12s}'*2+'\n')\
                .format('# STEP', 'RC', 'FEG')
        rcs = tot_rcs[:,i]
        thfos = tot_thfos[:,i]
        for k, (rc, thfo) in enumerate(zip(rcs, thfos)):
            content += ('{:<8d}'+'{:>12.4f}'*2+'\n')\
                .format(k+1, rc, thfo)

        datname = 'THFO-' + str(i+1) + '.dat'
        with open(datname, 'w') as writer:
            writer.write(content)

        print('Successfully write ThermoDynamic Forces to %s.' %datname)

    return 


def plot_lambs(tol=1000, paras={}):
    # setting
    tol = 1000 # withdraw 1000 steps
    paras = {'R': 0, 'A': 1} # indices and names of paras

    # read report
    tot_cons, tot_lamb = read('REPORT', tol)

    # first para
    fig, ax = plt.subplots(1, 1, figsize=(12,9))

    #cons = tot_cons[:,0]
    for name, i in paras.items():
        lamb = tot_lamb[:,i]
        ax.plot(range(len(lamb)), lamb, label=name)

    plt.legend()
    plt.savefig('sg.png')


def plot_mean_lamb(ncons, tol, para, intv=10):
    # setting
    nrelax = 50
    print('Drop first %d steps for relaxation.' %nrelax)
    
    # calc mean
    tot_cons, tot_lamb, tot_zdet, tot_zg = read('REPORT', ncons, tol)
    lambs = tot_lamb[nrelax:,para-1]
    zdets = tot_zdet[nrelax:,para-1]
    zgs = tot_zg[nrelax:,para-1]

    if ncons == 1:
        print('free energy gradient equals lambda')
    elif ncons > 1:
        print('free energy gradient equals <Z*(lambda+GkT)>/<Z-1/2>')
    else:
        raise ValueError('Incorrect number of constraints.')

    # write lambs to lamb.dat
    #lines = ['%.4f' %l for l in lamb]
    #with open('lamb.dat', 'w') as writer:
    #    writer.write('# Lambdas of Selected Constraint\n')
    #    writer.write('\n'.join(lines))

    intervals = np.arange(intv,len(lambs)+intv,intv)

    avgs, stds = [], []
    for i in intervals:
        dat = lambs[:i]
        avg = np.mean(zgs[:i]) / np.mean(zdets[:i]) # free energy gradient
        avgs.append(avg)
        stad = np.std(dat)
        stds.append(stad)

    # write data to BM-STAT.dat
    content = '# STEP AVG STD \n'
    for i, avg, std in zip(intervals, avgs, stds):
        content += '{:<8d}{:<8.4f}{:<8.4f}\n'. format(i, avg, std)
    with open('BM-STAT.dat', 'w') as writer:
        writer.write(content)

    # plot
    fig, ax = plt.subplots(1, 1, figsize=(12,9))
    ax.set_title('Free Energy Gradient', fontsize=20, fontweight='bold')

    ax.plot(intervals, [0]*len(avgs), color='k', linewidth=2)
    ax.plot(range(len(lambs)), lambs, color='r')

    ax.scatter(intervals, avgs, label='Average')
    ax.scatter(intervals, stds, label='Standard')

    plt.legend()
    plt.savefig('BM.png')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--ncons', type=int, \
            help='Number of Constraints')
    parser.add_argument('-nf', '--nframe', type=int, \
            help='Number of Steps to Average')

    args = parser.parse_args()

    read('TFILOG', args.ncons, args.nframe)
