#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os

import argparse

import numpy as np

import matplotlib as mpl
mpl.use('Agg') #silent mode
from matplotlib import pyplot as plt

MAXLINE = 10000

class SplineRepulsion():
    def __init__(self, npoints, cutoff, begrep, splrep, endrep):
        #
        self.npoints, self.cutoff = npoints, cutoff
        self.begrep = begrep
        self.splrep = splrep
        self.endrep = endrep

        #
        bounds = np.zeros(npoints+2)
        bounds[0], bounds[npoints+1] = 0., cutoff
        bounds[1:npoints] = splrep[:,0]
        bounds[npoints] = endrep[0]

        self.bounds = bounds

        return

    def calc_rep(self,r):
        """"""
        for i in range(self.npoints+1):
            b1, b2 = self.bounds[i], self.bounds[i+1]
            if b1 <= r < b2:
                if i == 0:
                    rep = self.calc_begrep(r,*self.begrep)
                elif 0 < i < self.npoints:
                    cur_spl = self.splrep[i-1,:]
                    coefs = np.zeros(5)
                    coefs[0] = cur_spl[0]
                    coefs[1:] = cur_spl[2:]
                    rep = self.calc_splrep(r,*coefs)
                elif i == self.npoints:
                    cur_spl = self.endrep
                    coefs = np.zeros(7)
                    coefs[0] = cur_spl[0]
                    coefs[1:] = cur_spl[2:]
                    rep = self.calc_endrep(r,*coefs)
                break
        else:
            raise ValueError('Distance %12.8f Not in Bound...' %r)

        return rep

    @staticmethod
    def calc_begrep(r,a1,a2,a3):
        rep = np.exp(-a1*r+a2) + a3

        return rep

    @staticmethod
    def calc_splrep(r,r0,c0,c1,c2,c3):
        rep = c0 + c1*(r-r0) + c2*(r-r0)**2 + c3*(r-r0)**3

        return rep

    @staticmethod
    def calc_endrep(r,r0,c0,c1,c2,c3,c4,c5):
        rep = c0 + c1*(r-r0) + c2*(r-r0)**2 + c3*(r-r0)**3 \
                + c4*(r-r0)**4 + c5*(r-r0)**5

        return rep


def plot_spline(skf='Pt-Pt.skf', pic='spl.png'):
    """Plot the Spline Repulsive Potential..."""
    # read spline data
    fopen = open(skf, 'r')
    for i in range(MAXLINE):
        line = fopen.readline()
        if line:
            if line.startswith('Spline'):
                # points and cutoff
                line = fopen.readline()
                data = line.strip().split()
                npoints, cutoff = int(data[0]), float(data[1])

                # exp rep
                # exp(-a1*r+a2)+a3
                line = fopen.readline()
                data = line.strip().split()
                begrep = np.array(data, dtype=float)

                # spline
                # c0+c1(r-r0)+c2(r-r0)**2+c3(r-r0)**3
                splrep = []
                for j in range(npoints-1):
                    line = fopen.readline()
                    data = line.strip().split()
                    splrep.append(data)
                splrep = np.array(splrep, dtype=float)

                # end
                line = fopen.readline()
                data = line.strip().split()
                endrep = np.array(data, dtype=float)

        else:
            # end of the file
            break
    fopen.close()

    # init spline
    sprp = SplineRepulsion(npoints,cutoff,begrep,splrep,endrep)

    # plot
    rs = np.linspace(0.,sprp.cutoff-0.01,1000)
    reps = []
    for r in rs:
        reps.append(sprp.calc_rep(r))
    skf_name = os.path.basename(skf).split('.')[0]

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12,8))
    ax.set_title(r'$%s$ Spline Repulsive Potential' %skf_name, \
            fontsize=24, fontweight='bold')
    ax.set_xlabel(r'$r$ / Bohr', fontsize=16)
    ax.set_ylabel(r'$V_{rep}(r)$ / Hartree', fontsize=16)

    ax.plot(rs, reps, color='k')

    plt.legend()

    plt.savefig(pic)

if __name__ == '__main__':
    """
    parser = argparse.ArgumentParser()

    parser.add_argument('-f', '--skf', help='Slater-Koster File')
    parser.add_argument('-p', '--pic', \
            default='spl.png', help='Spline Repulsive Potential Figure')

    args = parser.parse_args()
    """

    #plot_spline(args.skf, args.pic)
    plot_spline()
