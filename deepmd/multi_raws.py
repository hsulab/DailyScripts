#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os 
import subprocess
import argparse 

import numpy as np 

from joblib import Parallel, delayed

RAW_TO_SET='/users/40247882/projects/ucdp/copper/utils/raw_to_set.sh'

def run_command(command, cwd, timeout=120):
    proc = subprocess.Popen(
        command, shell=True, cwd=cwd, 
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, 
        encoding = 'utf-8'
    )
    errorcode = proc.wait(timeout=timeout) # 10 seconds 
    if errorcode:
        print(cwd)
        raise ValueError('Error in command %s.' %(command))

    return proc

if __name__ == '__main__':
    """"""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-s', '--symbol', 
        default='Cu', help='chemical symbol'
    )
    parser.add_argument(
        '-ns', '--nsets', type=int, 
        default=5, help='chemical symbol'
    )
    parser.add_argument(
        '-nj', '--njobs', type=int,
        default=4, help='number of processes'
    )

    args = parser.parse_args()

    raw_dirs = []
    for dname in os.listdir('./'):
        if dname.startswith(args.symbol):
            raw_dirs.append(dname)
    #raw_dirs.sort()

    with open('raw.log','w') as fopen:
        fopen.write('')

    commands = []
    for dname in raw_dirs:
        content = 'System %s\n' %dname
        # number of frames 
        proc = run_command('cat box.raw | wc -l', dname)
        nframes = int(proc.stdout.readlines()[0].strip())
        n_per_set = np.ceil(nframes/args.nsets)
        #if n_per_set == 0:
        #    n_per_set += 1
        n_per_set = int(n_per_set)
        content += 'nframes %d\n' %nframes
        content += 'nframe per set %d\n' %n_per_set
        ## raw to set
        cur_mand = RAW_TO_SET + ' %d' %n_per_set
        commands.append(cur_mand)
        proc = run_command(cur_mand, dname)
        content += ''.join(proc.stdout.readlines())
        content += '\n\n' 
        print(content)
        with open('raw.log','a') as fopen:
            fopen.write(content)

    #print('now actual run...')
    #print('using num of jobs: ', args.njobs)
    #Parallel(n_jobs=args.njobs)(delayed(run_command)(cur_mand, dname) for dname, cur_mand in zip(raw_dirs, commands))

    pass
