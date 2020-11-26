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

def read_report(filename='REPORT', ncons=2, tolstep=3):
    # get file object
    reader = open(filename, 'r')
        
    nstep, flag = 0, True
    tot_cons, tot_lamb = [], []
    cur_cons, cur_lamb = [], []
    tot_zdet, tot_zg = [], []
    cur_zdet, cur_zg = [], []
    while flag:
        # find MD step
        line = reader.readline().strip()
        if line.startswith('MD step No.'):
            # count and reset data
            nstep += 1
            cur_cons, cur_lamb = [], []
            cur_zdet, cur_zg = [], []
            irbm = False # If Read Blue Moon

        # read section
        if line.startswith('>Const_coord'):
            for i in range(ncons):
                line = reader.readline().strip()
                val = float(line.split()[2])
                cur_cons.append(val)
            tot_cons.append(cur_cons)

        if line.startswith('>Blue_moon'):
            line = reader.readline().strip() # lambda |z|^(-1/2) GkT ... 
            for i in range(ncons):
                line = reader.readline().strip()
                vals = [float(val) for val in line.split()[1:]]
                cur_lamb.append(vals[0])
                cur_zdet.append(vals[1])
                cur_zg.append(vals[3])
            tot_lamb.append(cur_lamb)
            tot_zdet.append(cur_zdet)
            tot_zg.append(cur_zg)
            irbm = True

        # check step
        if nstep == tolstep and irbm == True:
            flag = False

    reader.close()
    print('Successfully read CONSTRAINTS and LAMBDAS in REPORT.')

    # check output
    tot_cons = np.array(tot_cons)
    tot_lamb = np.array(tot_lamb)
    tot_zdet = np.array(tot_zdet)
    tot_zg = np.array(tot_zg)

    # write constraints and lambdas to BM.dat
    for i in range(ncons):
        content = ('{:<8s}'+'{:<12s}'*5+'\n')\
                .format('# STEP', 'CV', 'lamb', '|z|^(-1/2)', '||*(lamb+G)', 'Gradient')
        cons = tot_cons[:,i]
        lambs = tot_lamb[:,i]
        zdets = tot_zdet[:,i]
        zgs = tot_zg[:,i]
        for k, (con, lamb, zdet, zg) in enumerate(zip(cons, lambs, zdets, zgs)):
            content += ('{:<8d}'+'{:>12.4f}'*5+'\n')\
                .format(k+1, con, lamb, zdet, zg, zg/zdet)

        datname = 'BM-' + str(i+1) + '.dat'
        with open(datname, 'w') as writer:
            writer.write(content)

        print('Successfully write LAMBDAS to %s.' %datname)

    return tot_cons, tot_lamb, tot_zdet, tot_zg


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
    plt.savefig('BM-'+str(para)+'.png')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--ncons', type=int, \
            help='Number of Constraints')
    parser.add_argument('-nf', '--nframe', type=int, \
            help='Number of Steps to Average')
    parser.add_argument('-p', '--parameter', type=int, \
            help='Constrained Parameter Index')

    args = parser.parse_args()

    read_report('REPORT', args.ncons, args.nframe)
