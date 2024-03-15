#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse

import numpy as np

import matplotlib
matplotlib.use('Agg') #silent mode
import matplotlib.pyplot as plt

from pymatgen.io.vasp import Locpot

def e_vacuum():
    slab_locpot = Locpot.from_file("LOCPOT")
    with open('POSCAR','r') as f:
        content = f.readlines()
    z_length = eval(content[4].split()[-1])
    slab_zlp = slab_locpot.get_average_along_axis(2)
    print("LOCPOT z-axis data: ", slab_zlp.shape)
    vac_range = int(3/z_length*len(slab_zlp))
    print(vac_range)
    #return np.mean(slab_zlp[150:200])
    return np.mean(slab_zlp[175])

def get_fermi():
    with open ('OUTCAR','r') as f:
        context = f.readlines()
    target_line = []
    for l in context:
        if l.find('E-fermi') == -1:
            continue
        else:
            target_line.append(l)
    fermi = target_line[-1].split()[2]
    # print('Fermi Energy {fermi}')
    return eval(fermi)

def plot_zpotential():
    slab_locpot = Locpot.from_file("LOCPOT")
    slab_zlp = slab_locpot.get_average_along_axis(2)
    plt.plot(slab_zlp)
    plt.savefig('zpotential.png')

# Main Function
plot_zpotential()
vacuum_energy = e_vacuum()
fermi_energy = get_fermi()
content = "===== Workfunction Overview =====\n"
content += "Vacuum Energy: %.4f\n" %vacuum_energy
content += "Ferimi Energy: %.4f\n" %fermi_energy
content += "Workfunction : %.4f\n" %(vacuum_energy-fermi_energy)
print(content)

if __name__ == '__main__':
    pass
