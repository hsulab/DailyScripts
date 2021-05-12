#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os

import argparse

import numpy as np
import pandas as pd

#import pyblock as pb

import matplotlib as mpl
mpl.use('Agg') #silent mode
from matplotlib import pyplot as plt

"""
Author: Jiayan Xu
"""

NLINTRO = 50 # number of lines in intro part
NLSTEP = 200
MAXCON = 50

MAXFRAME = 100000

BINWIDTH = 0.005

class Report(object):
    def __init__(self, fname='REPORT'):
        # self.ftype = 'US' # filetype US or BM
        #self.fname = 'REPORT-US' # REPORT
        if os.path.exists(fname):
            self.fname = fname
        else:
            raise ValueError('%s does not exist.' %fname)

    def parse_report(self, datfile, nframe, ndrop):
        reader = open(self.fname, 'r')

        reader, ftype, ncon = self.parse_intro(reader)
        reader, ncon = self.parse_fstep(reader, ftype)

        if ftype == 'US':
            #if os.path.exists(datfile):
            #    print('%s already exists.' %datfile)
            #    return
            # read reaction coordinates
            tot_coords = []
            for i in range(nframe):
                reader, coords, tend = self.parse_us_step(reader, ncon)
                if tend:
                    print('%s only has %d steps.' %(self.fname, i))
                    break
                tot_coords.append(coords)
            # print(tot_coords)
            tot_coords = np.array(tot_coords)

            # write to full/drop dat
            cmins, cmaxs = [], []
            content = '# step coord\n'
            for n in range(ncon):
                cmin, cmax = np.min(tot_coords[:,n]), np.max(tot_coords[:,n])
                content += '# %d min: %12.8f max: %12.8f\n' \
                        %(n+1, cmin, cmax)
                cmins.append(cmin)
                cmaxs.append(cmax)
            self.cmins, self.cmaxs = cmins, cmaxs

            con_list = []
            for i, coords in enumerate(tot_coords):
                con_list.append(('{:<8d}'+'{:>12.8f}'*ncon)\
                        .format(i+1, *coords))
            if datfile:
                bname = (datfile.split('.'))[0] # /root/home/US.dat
            else:
                bname = 'US'
            with open(bname+'.dat', 'w') as writer:
                writer.write(content+'\n'.join(con_list))
            if ndrop > 0:
                with open(bname+'-drop.dat', 'w') as writer:
                    writer.write(content+'\n'.join(con_list[ndrop:]))
            data = tot_coords
        elif ftype == 'BM':
            tot_coords, tot_lambs = [], []
            tot_zdets, tot_gkts, tot_zgs = [], [], []
            for i in range(nframe):
                reader, coords, lambs, zdets, gkts, zgs, tend = \
                        self.parse_bm_step(reader, ncon)
                if tend:
                    print('%s only has %d steps.' %(self.fname, i))
                    break
                tot_coords.append(coords)
                tot_lambs.append(lambs)
                tot_zdets.append(zdets)
                tot_gkts.append(gkts)
                tot_zgs.append(zgs)
            tot_coords = np.array(tot_coords)
            tot_lambs = np.array(tot_lambs)
            tot_zdets = np.array(tot_zdets)
            tot_gkts = np.array(tot_gkts)
            tot_zgs = np.array(tot_zgs)

            # write to file
            if datfile:
                bname = (datfile.split('.'))[0] # /root/home/US.dat
            else:
                bname = 'BM-'

            for n in range(ncon):
                content = ('# '+'{:<8s}'*6+'\n')\
                    .format('Step', 'Coord', 'Lamb', '|Z|', 'GkT', '|Z|^*()')
                for i, (coord, lamb, zdet, gkt, zg) in enumerate(
                    zip(tot_coords[:,n],tot_lambs[:,n],\
                        tot_zdets[:,n],tot_gkts[:,n],tot_zgs[:,n])):
                    content += ('{:<8d}'+'  {:>8.4f}'*5+'\n')\
                        .format(i+1, coord, lamb, zdet, gkt, zg)
                datfile = bname+str(n+1)+'.dat'
                with open(datfile, 'w') as writer:
                    writer.write(content)

            data = [tot_coords, tot_lambs, tot_zdets, tot_gkts, tot_zgs]

        reader.close()
        print('Read %s and Write %s.' %(self.fname, datfile))

        return ftype, data

    @staticmethod
    def parse_intro(reader):
        """parse intro lines in REPORT"""
        for i in range(NLINTRO):
            line = reader.readline()
            if line.startswith('                MDALGO'):
                tag = (line.strip().split())[-1]
                if tag.isdigit():
                    tag = int(tag)
                elif tag == '**':
                    tag = 21
                else:
                    raise ValueError('MDALGO should be an integer.')

                if tag == 2 or tag == 26:
                    ftype = 'BM'
                elif tag == 21 or tag == 27:
                    ftype = 'US'
                else:
                    raise ValueError('Unsupported MDALGO REPROT.')
            #print(line, end='')
            if line.startswith('   original number of atomic DOF'):
                ndof_org = int((line.strip().split())[-1])
                #
                line = reader.readline()
                ncon = int((line.strip().split())[-1])
                line = reader.readline()
                ndof_act = int((line.strip().split())[-1])
                #print(ndof_org, ncon, ndof_act)
                #print(ftype)
                break
        else:
            raise ValueError('Intro lines in REPORT may be incorrect.')

        return reader, ftype, ncon

    @staticmethod
    def parse_fstep(reader, ftype):
        curln = reader.tell()
        """read first step"""
        if ftype == 'US':
            for i in range(NLSTEP):
                line = reader.readline()
                if line.startswith('  >Metadynamics'):
                    coords = []
                    for j in range(MAXCON):
                        line = reader.readline()
                        if line.startswith('   fic_p>'):
                            coord = float((line.strip().split())[-1])
                        else:
                            ncon = j
                            break
                    else:
                        raise ValueError('Too many lines in >Metadynamics.')
                if line.startswith('           RANDOM_SEED'):
                    break
            else:
                raise ValueError('Too many lines in =MD STEP=.')
        elif ftype == 'BM':
            for i in range(NLSTEP):
                line = reader.readline()
                if line.startswith('  >Const_coord'):
                    coords = []
                    for j in range(MAXCON):
                        line = reader.readline()
                        if line.startswith('   cc>'):
                            # cc> tag val val tol
                            coord = float((line.strip().split())[2])
                        else:
                            ncon = j
                            break
                    else:
                        raise ValueError('Too many lines in >Metadynamics.')
                if line.startswith('           RANDOM_SEED'):
                    break
            else:
                raise ValueError('Too many lines in =MD STEP=.')

        reader.seek(curln)

        return reader, ncon

    @staticmethod
    def parse_us_step(reader, ncon):
        """read US step"""
        tend = False # if read to the end line
        coords = []
        for i in range(NLSTEP):
            line = reader.readline()
            if not line:
                tend = True
                break
                # raise ValueError('dsad')
            # print(line, end='')
            if line.startswith('  >Metadynamics'):
                for j in range(ncon):
                    line = reader.readline()
                    coord = float((line.strip().split())[-1])
                    coords.append(coord)
            #if line.startswith('           RANDOM_SEED'):
            if line.startswith('           RANDOM_SEED'):
                # line = reader.readline()
                break
        else:
            raise ValueError('Too many lines in =MD STEP=.')

        return reader, coords, tend

    @staticmethod
    def parse_bm_step(reader,ncon):
        """read BM step"""
        tend = False # if read to the end line
        coords, lambs, zdets, gkts, zgs = [], [], [], [], []

        for i in range(NLSTEP):
            line = reader.readline()
            if not line:
                tend = True
                break
                # raise ValueError('dsad')
            # print(line, end='')
            if line.startswith('  >Const_coord'):
                for j in range(ncon):
                    line = reader.readline()
                    coord = float((line.strip().split())[2])
                    coords.append(coord)
            if line.startswith('  >Blue_moon'):
                line = reader.readline() # comment line
                for j in range(ncon):
                    line = reader.readline()
                    vals = [float(val) for val in (line.strip().split())[1:]]
                    lambs.append(vals[0])
                    zdets.append(vals[1])
                    gkts.append(vals[2])
                    zgs.append(vals[3])
            if line.startswith('           RANDOM_SEED'):
                # line = reader.readline()
                break
        else:
            raise ValueError('Too many lines in =MD STEP=.')

        return reader, coords, lambs, zdets, gkts, zgs, tend

def read_input(infile, ndrop):
    with open(infile, 'r') as reader:
        lines = reader.readlines()
    lines = [line.strip().split() \
            for line in lines if not line.strip().startswith('#')]
    datfiles = [line[0] for line in lines]

    rcmins, rcmaxs = [], []
    for datfile in datfiles:
        usdir = os.path.dirname(datfile)
        print(usdir)
        report = Report(os.path.join(usdir,'REPORT'))
        report.parse_report(datfile, nframe=MAXFRAME, ndrop=ndrop)
        rcmins.append(report.cmins)
        rcmaxs.append(report.cmaxs)
    rcmins, rcmaxs = np.array(rcmins), np.array(rcmaxs)

    ncon = rcmins.shape[1]
    for i in range(ncon):
        rcmin = np.min(rcmins[:,i])
        rcmax = np.max(rcmaxs[:,i])
        print('Reaction Coordinate %4d Range %.4f to %.4f' \
                %(i+1, rcmin, rcmax))

    return

def read_input2(infile, ndrop):
    with open(infile, 'r') as reader:
        lines = reader.readlines()
    lines = [line.strip().split() \
            for line in lines if not line.strip().startswith('#')]
    datfiles = [line[0] for line in lines]

    rcmins, rcmaxs = [], []
    for datfile in datfiles:
        usdir = os.path.dirname(datfile)
        report = Report(os.path.join(usdir,'REPORT'))
        report.parse_report(datfile, nframe=MAXFRAME, ndrop=ndrop)
        rcmins.append(report.cmins)
        rcmaxs.append(report.cmaxs)
    rcmins, rcmaxs = np.array(rcmins), np.array(rcmaxs)

    ncon = rcmins.shape[1]
    for i in range(ncon):
        rcmin = np.min(rcmins[:,i])
        rcmax = np.max(rcmaxs[:,i])
        print('Reaction Coordinate %4d Range %.4f to %.4f' \
                %(i+1, rcmin, rcmax))

    return

def read_usdat(datfile):
    """"""
    with open(datfile, 'r') as reader:
        lines = reader.readlines()
    lines = [line.strip().split() for line in lines \
            if not line.startswith('#')]
    data = np.array(lines,dtype=float)

    return data

def plot_usdata(datfile, col):
    """"""
    data = read_usdat(datfile)

    step = data[:,0]

    for c in col:
        cval = data[:,c]

        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12,8))

        ax.set_title('Umbrella Sampling', fontsize=20, fontweight='bold')
        ax.set_xlabel('Step (Number)', fontsize=16)
        ax.set_ylabel('Collective Variable', fontsize=16)

        ax.plot(step, cval, color='b')

        nbins = int(np.ceil((np.max(cval)-np.min(cval))/BINWIDTH)+1)

        n, bins, p = ax.hist(cval, bins=nbins, \
                color='g', alpha=0.2, orientation='horizontal')
        centres = []
        for i, b in enumerate(bins[:-1]):
            centres.append((b+bins[i+1])/2.0)

        ax.plot(n, centres, color='r')

        plt.savefig('uscv-'+str(c)+'.png')

    return

def plot_bmdata(data, cols):
    lambs = data[1] # lambdas 
    dim = lambs.shape

    for col in cols:
        if col > dim[1]:
            raise ValueError('Not so many constraints.')
        data = lambs[:,col-1].ravel()

        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12,8))

        ax.set_title('Constrained Molecular Dynamics', \
                fontsize=20, fontweight='bold')
        ax.set_xlabel('Step (Number)', fontsize=16)
        ax.set_ylabel(r'$\lambda$', fontsize=16)

        ax.plot(np.arange(len(data))+1, data, color='b')
        ax.plot(np.arange(len(data))+1, [0]*len(data), color='k', linestyle='--')

        plt.savefig('BM-'+str(col)+'.png')

    return

#def block_average(data, cols, ndrop):
#    lambs = data[1]
#    dim = lambs.shape
#
#    for col in cols:
#        if col > dim[1]:
#            raise ValueError('Not so many constraints.')
#        data = pd.Series(lambs[ndrop:,col-1].ravel())
#
#        (data_length, reblock_data, covariance) = pb.pd_utils.reblock(data)
#        block_info = reblock_data
#
#        data_sets = block_info.columns.get_level_values(0).unique()
#        for i, ci in enumerate(data_sets):
#            opt = block_info[block_info[(ci,'optimal block')] != ''].index.values
#            if opt:
#                opt_mean = block_info.loc[opt[0],(ci,'mean')]
#                print('BM-%d %5d/%5d MEAN %.4f' %(col,ndrop,len(data),opt_mean))
#                break
#        else:
#            print('No opt...')
#
#    return 


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--version', \
            action='version', version='%(prog)s 1.4')

    group = parser.add_mutually_exclusive_group()
    group.add_argument('-mt', '--meta', \
            action='store', help='METADATA file')
    group.add_argument('-rp', '--repo', nargs='?',\
            default=None, action='store', help='REPORT file')

    parser.add_argument('-nd', '--ndrop', type=int, nargs='?',\
            default=-1, help='Number of Steps Dropped')

    #parser.add_argument('-n', '--ncons', type=int, \
    #        help='Number of Constraints')
    #parser.add_argument('-nf', '--nframe', type=int, \
    #        help='Number of Steps to Average')
    parser.add_argument('-p', '--para', type=int, nargs='*',\
            default=[-1], help='plot the p-th collective variable')

    args = parser.parse_args()

    #read_report('REPORT', args.ncons, args.nframe)
    #report = Report('REPORT-US')
    #report.parse_report(nframe=3)

    #read_input('METADATA', ndrop=1000)
    if args.repo:
        report = Report(args.repo)
        ftype, data = report.parse_report(None, nframe=MAXFRAME, ndrop=args.ndrop)
        if args.para[0] != -1:
            if ftype == 'US':
                plot_usdata('US.dat', args.para)
            elif ftype == 'BM':
                plot_bmdata(data, args.para)
                if args.ndrop != -1:
                    # block_average(data, args.para, args.ndrop)
                    raise NotImplementedError('block average not supported')
            print('plot the p-th collective variable.')
    elif args.meta:
        print('Read METADATA ...')
        read_input(args.meta, args.ndrop)
    else:
        parser.print_help()
