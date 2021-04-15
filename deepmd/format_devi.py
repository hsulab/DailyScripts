#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np 

with open('model_devi.out') as fopen:
    lines = fopen.readlines() 

nsteps = len(lines) - 1 # the first line is comment 
print('number of steps %d' %nsteps)

devi_stat = [] # max_devi_e, min_devi_e, avg_devi_e, ...
atomic_devi = []
for i in range(1,nsteps+1):
    data = lines[i].strip().split() 
    print(len(data))
    devi_stat.append(data[1:7])
    atomic_devi.append(data[7:])

assert len(devi_stat) == len(atomic_devi) 

devi_stat = np.array(devi_stat, dtype=np.float32)
print(devi_stat.shape)
atomic_devi = np.array(atomic_devi, dtype=np.float32)
print(atomic_devi.shape)
atomic_devi = np.reshape(atomic_devi, (nsteps,-1))
print(atomic_devi)


if __name__ == '__main__':
    pass
