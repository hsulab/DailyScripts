#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os

import argparse

import numpy as np

import matplotlib as mpl
mpl.use('Agg') #silent mode
from matplotlib import pyplot as plt

from mpl_toolkits.axes_grid1.inset_locator import inset_axes, InsetPosition, zoomed_inset_axes, mark_inset

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

def read_spline(skf):
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

    return sprp


def plot_spline(skf='Pt-Pt.skf', skf2=None, rmin=1.0, pic='spl.png'):
    """Plot the Spline Repulsive Potential..."""
    # read spline, turn into spline object
    SP_rep1 = read_spline(skf)

    # generate data
    rs = np.linspace(0.,SP_rep1.cutoff-0.01,1000)
    reps = []
    for r in rs:
        reps.append(SP_rep1.calc_rep(r))
    skf_name = os.path.basename(skf).split('.')[0]

    rs = np.array(rs)
    reps = np.array(reps)

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12,8))
    ax.set_title(r'$%s$ Spline Repulsive Potential' %skf_name, \
            fontsize=24, fontweight='bold')
    ax.set_xlabel(r'$r$ / Bohr', fontsize=20)
    ax.set_ylabel(r'$V_{rep}(r)$ / Hartree', fontsize=20)

    skf1_curve, = ax.plot(rs, reps, \
            color='g', linestyle='-', linewidth=2., \
            label='Skf-1')

    # inset figure
    ax2 = plt.axes([0,0,1,1])
    ip = InsetPosition(ax, [0.4,0.2,0.5,0.5])
    ax2.set_axes_locator(ip)
    mark_inset(ax, ax2, loc1=1, loc2=3, fc="none", ec='0.5')

    #ax2 = zoomed_inset_axes(ax, 1, loc=4)
    r_min, r_max = rmin, SP_rep1.cutoff
    indices = np.where((rs>r_min) & (rs<r_max))
    ax2.plot(rs[indices], reps[indices], color='g', linestyle='-', linewidth=2.)

    # skf2 for comparision
    if skf2:
        SP_rep2 = read_spline(skf2)

        # generate data
        rs = np.linspace(0.,SP_rep2.cutoff-0.01,1000)
        reps = []
        for r in rs:
            reps.append(SP_rep2.calc_rep(r))

        rs = np.array(rs)
        reps = np.array(reps)

        skf2_curve, = ax.plot(rs, reps, \
                color='orange', linestyle='--', linewidth=2., \
                label='Skf-2')
        ax2.plot(rs[indices], reps[indices], color='orange', linestyle='--', linewidth=2.)
        plt.legend(handles=[skf1_curve,skf2_curve])
    else:
        plt.legend(handles=[skf1_curve,])

    plt.savefig(pic)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-f', '--skf', required=True,\
            help='Slater-Koster File')
    parser.add_argument('-f2', '--skf2', default=None, \
            help='Second Slater-Koster File for Comparision')
    parser.add_argument('-p', '--pic', \
            default='spl.png', help='Spline Repulsive Potential Figure')
    parser.add_argument('-rmin', '--radius_min', type=float,\
            default=1.0, help='Minimum Radius for Zoom')

    args = parser.parse_args()

    plot_spline(args.skf, args.skf2, args.radius_min, args.pic)
    #plot_spline()
