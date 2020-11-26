#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os

import argparse

import numpy as np

import matplotlib as mpl
mpl.use('Agg') #silent mode
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

ELEMENT_COLOR = {'default': '#000000', 'H': '#FFFFF0', 'C': '#A9A9A9', 
'O': '#FF0000', 'S': '#FEC615', 'Cd': '#645403', 'Pt': '#0000CD'}
ORBITAL_COLOR = {'s': '#4b0101', 'p': '#01366a', 'd': '#0a461e', 'f': '#751973'}

ORBITALS = ['s', 'py', 'pz', 'px', \
    'dxy', 'dyz', 'dz2', 'dxz', 'dx2', \
    'f-3', 'f-2', 'f-1', 'f0', 'f1', 'f2', 'f3']

class AtomDos(object):
    """
    Description:
        A class for a single atom in DOSCAR.
    ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
    element attr.  [str]   the element symbol
    number  attr.  [int]   the index, start from 1
    _nedos  attr.  [int]   number of bands
    _doses  attr.  [dict]  {'up': [], 'dw': []}
    ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
    parse_block  meth.
    get_dos      meth.
    ===== ===== ===== ===== ===== ===== ===== ===== ===== =====

    """
    def __init__(self, element, number, block):
        self._element = element
        self._number = number
        self._nedos, self._doses = self.parse_block(block) # dict

    @property
    def element(self):
        return self._element

    @property
    def number(self):
        return self._number

    def parse_block(self, block):
        # energy s py pz px dxy dyz dz2 dxz dx2 f-3 f3
        #print(block[0])
        emax, emin, nedos, efermi, scale = block[0] # these are all strings
        nedos = int(nedos)

        dos_data = np.array(block[1:], dtype=float)

        doses = dos_data[:,1:]

        # s 1 sp 4 spd 9 spdf 16
        # check orbitals
        norbitals = len(doses[0])

        if norbitals == 1 or norbitals == 4 or norbitals == 9 or norbitals == 16:
            orbitals = ORBITALS[:norbitals]
            doses_dict = {}
            for i, orbital in enumerate(orbitals):
                doses_dict[orbital] = doses[:,i]
        elif norbitals == 2 or norbitals == 8 or norbitals == 18 or norbitals == 32:
            # spin-polarized
            orbitals = ORBITALS[:int(norbitals/2)]
            doses_up, doses_dw = [], []
            for i in range(norbitals):
                if i%2 == 0:
                    doses_up.append(doses[:,i])
                else:
                    doses_dw.append(doses[:,i])

            #doses_dict_up, doses_dict_dw = {}, {}
            doses_dict = {}
            for i, orbital in enumerate(orbitals):
                doses_dict[orbital+'_up'] = doses_up[i]
                doses_dict[orbital+'_dw'] = doses_dw[i]

        return nedos, doses_dict

    def get_dos(self, orbital=None, spin=None):
        """
        >>> AtomDos.get_dos('s') # get s_up-s_dw orbital
        >>> AtomDos.get_dos('p', 'up') # get s_up orbital

        """
        doses = np.zeros(self._nedos)
        if orbital in ['s', 'p', 'd', 'f']:
            if spin == 'up' or spin == 'dw':
                spins = [spin]
            else:
                spins = ['up', 'dw']
            for o, d in self._doses.items():
                if o.startswith(orbital):
                    if spin:
                        if spin in o:
                            doses += np.fabs(d)
                    else:
                        doses += np.fabs(d)
                else:
                    pass
        elif orbital in ORBITALS:
            if spin == 'up' or spin == 'dw':
                orbital += '_' + spin
            for o, d in self._doses.items():
                if orbital in o:
                    doses += d
        else:
            raise ValueError('Unknown Orbital Type.')

        return doses

    def __repr__(self):
        return 'DOS %s%s (%d)' %(self.element, self.number, len(self._doses))

class DosCar(object):
    """
    Description:
        A class for a single atom in DOSCAR.
    ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
    _doses
    ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
    ===== ===== ===== ===== ===== ===== ===== ===== ===== =====

    """
    def __init__(self, filename='DOSCAR'):
        # doscar name
        self.filename = filename

        # get atom list
        with open('POSCAR', 'r') as reader:
            lines = reader.readlines()
        elements = lines[5].strip().split()
        numbers = [int(i) for i in lines[6].strip().split()]

        atoms = []
        for e, n in zip(elements, numbers):
            atoms.extend([e]*n)

        self._atoms = atoms

        # read DOSCAR
        self._doses = self.readfile()

    def readfile(self):
        with open(self.filename, 'r') as reader:
            lines = reader.readlines()
        lines = [line.strip().split() for line in lines]

        # parse data
        # general info
        natoms = int(lines[0][0])
        a_omega, a_norm, b_norm, c_norm, potim = np.array(lines[1], dtype=float)
        system_name = lines[4]

        # total dos
        emax, emin, nedos, e_fermi, scale = np.array(lines[5], dtype=float)
        nedos = int(nedos)
        if len(lines[6]) == 3:
            ISPIN = 1
        elif len(lines[6]) == 5:
            ISPIN = 2
        else:
            raise ValueError('Unrecognized DOSCAR Format.')

        self._nedos = int(nedos)

        dos_data = np.array(lines[6:6+nedos], dtype=float) # including energies

        self._energies = np.array(dos_data[:,0], dtype=float)

        # total dos, integrated total dos
        if ISPIN == 1:
            self.tot_dos = np.array(dos_data[:,1], dtype=float)
            self.tot_idos = np.array(dos_data[:,2], dtype=float)
        elif ISPIN == 2:
            self.tot_dos = np.array(dos_data[:,1:3], dtype=float)
            self.tot_idos = np.array(dos_data[:,3:], dtype=float)

        # projected dos
        atoms_dos = []
        for i, atom in enumerate(self._atoms):
            start, end = 6+nedos+i*(nedos+1), 6+nedos+(i+1)*(nedos+1)
            dos_data = lines[start:end]
            atoms_dos.append(AtomDos(atom, i+1, dos_data))

        return atoms_dos

    def get_energy(self):
        """band energy"""
        return self._energies

    def get_tdos(self):
        return self.tot_dos

    def get_pdos(self, element, number=None, orbital=None, spin=None):
        """
        # all S atoms and all orbitals
        >>> DosCar.get_pdos('S') 
        # selected S atoms and all orbitals
        >>> DosCar.get_pdos('S', [1]) 
        # all S atoms and selected orbitals
        >>> DosCar.get_pdos('S', [], ['s']) 
        # all S atoms and all orbitals and selected spin
        >>> DosCar.get_pdos('S', [], [], ['up']) ***
        """
        # sum pdos 
        pdos = np.zeros(self._nedos)
        for ados in self._doses:
            if ados.element == element:
                for s in spins:
                    pdos += ados.get_dos(orbital, s)

        return pdos

    def get_element_pdos(self, element, orbital=None, spin=None):
        """
        """
        pdos = np.zeros(self._nedos)
        for ados in self._doses:
            if ados.element == element:
                pdos += ados.get_dos(orbital, spin) # orbital spin

        return pdos


def plot_dos(dos_data_file, dos_fig_format):
    #
    fig_name = 'DOS.png'
    print('Reas DOS data of %s. Write figure to %s.')

    # read fermi energy
    with open('DOSCAR', 'r') as reader:
        lines = reader.readlines()
    e_fermi = float(lines[5].strip().split()[3])

    # read DOSCAR
    doscar = DosCar()
    energies = doscar.get_energy()

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8,6))
    #ax.set_title('Projected Density of States', \
    #        fontsize=20, fontweight='bold')

    ax.set_xlabel('Energy / eV', fontsize=16)
    ax.xaxis.set_major_locator(MultipleLocator(5.0))
    ax.xaxis.set_minor_locator(MultipleLocator(1.0))

    ax.set_ylabel('Density of States', fontsize=16)

    ax.plot(energies, len(energies)*[0], color='k')

    doses = doscar.get_element_pdos('Cd', 'd', 'up')
    ax.plot(energies, doses, color='g', label=r'Cd $d$')

    doses = doscar.get_element_pdos('Cd', 's', 'up')
    ax.plot(energies, doses, color='r', label=r'Cd $s$')
    
    doses = doscar.get_element_pdos('S', 'p', 'up')
    ax.plot(energies, doses, color='y', label=r'S $p$')

    doses = doscar.get_tdos()[:,0]
    ax.plot(energies, doses, color='k', label=r'Total DOS')

    plt.legend()

    plt.savefig(fig_name, format='png')

if __name__ == '__main__':
    #parser = argparse.ArgumentParser()

    #parser.add_argument('-f', '--datafile', \
    #        default='PDOS_USER.dat', help='PDOS Data File')
    #parser.add_argument('-o', '--out', \
    #        default='png', help='DOS Figure Format')
    #parser.add_argument('-d', '--delete', \
    #        help='Delete Element\'s DOS')

    #args = parser.parse_args()

    #plot_dos(args.datafile, args.out)
    plot_dos('','')
    pass
