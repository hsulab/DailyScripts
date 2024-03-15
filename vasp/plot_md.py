#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re

import numpy as np

import matplotlib
matplotlib.use('Agg') #silent mode
import matplotlib.pyplot as plt
plt.style.use("presentation")

#from ase.io import read, write

def parse_oszicar(oszicar):
    # read lines
    with open(oszicar, "r") as fopen:
        lines = fopen.readlines()

    # match lines
    data_lines = []
    for line in lines:
        if re.match("^[\t ]*[0-9]* T=*", line):
            data_lines.append(line)

    # index T E  F E0 EK SP SK
    data = []
    for line in data_lines:
        cur_data = [x for i, x in enumerate(line.strip().split()) if i%2 == 0]
        data.append(cur_data)
    data = np.array(data, dtype=float)

    return data

data = parse_oszicar("./OSZICAR")

steps = data[:, 0]
temperatures = data[:, 1]
energies = data[:, 4] # potential energy

#frames = read("./OUTCAR", ":")
#energies = [a.get_potential_energy() for a in frames]

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(16,12))
ax.set_title(
    "Molecular Dynamics"
)

l1, = ax.plot(steps, temperatures, c="r", label="temperature")
ax.set_xlabel("MD Step")
ax.set_ylabel("Temperature [K]")

ax2 = ax.twinx()
l2, = ax2.plot(steps, energies, label="potential energy")
ax2.set_ylabel("Energy [eV]")

plt.legend(handles=[l1, l2])

plt.tight_layout()
plt.savefig('md.png')

if __name__ == '__main__':
    pass
