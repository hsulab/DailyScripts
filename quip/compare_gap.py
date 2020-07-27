#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import ase.io
from ase import Atom, Atoms

from ase.optimize import BFGS
from ase.constraints import StrainFilter # optimize the unit cell while keeping scaled positiions fixed
from ase.constraints import FixAtoms

from quippy.potential import Potential
from quippy import descriptors

if __name__ == '__main__':
    # load GAP and frames
    gap = Potential(param_filename='./gap-4/GAP.xml')

    xyzPath = '/home/e89/e89/jyxu/projects/dftb-ml/std-dft/structures/reax_structures.xyz'
    reaction_frames = ase.io.read(xyzPath, ':') # ase atom frames

    structureList = [
            'O2', 'CO', 'CO2', # molecules
            'Surf', 'O_FCC', 'CO_Top', # adsorption
            'IS', 'TS_HCP', 'FS', # reaction
            ]

    # first, calculate standard dft energy 
    dft_energies = {} # using electronic free energy, in consistent with GAP
    gap_energies = {} # using electronic free energy, in consistent with GAP
    for i, sname in enumerate(structureList):
        atoms = reaction_frames[i]
        # DFT
        dft_energies[sname] = atoms.get_potential_energy(force_consistent=True)
        # GAP
        atoms.set_calculator(gap)
        gap_energies[sname] = atoms.get_potential_energy(force_consistent=True)

    # CO2 Formation Energy
    e_co2form = dft_energies['CO2'] - dft_energies['CO'] - 0.5*dft_energies['O2']
    e_oadfcc = dft_energies['O_FCC'] - dft_energies['Surf'] - 0.5*dft_energies['O2']
    e_coadtop = dft_energies['CO_Top'] - dft_energies['Surf'] - dft_energies['CO']
    e_fwbarr = dft_energies['TS_HCP'] - dft_energies['IS']
    e_bkbarr = dft_energies['TS_HCP'] - dft_energies['FS']
    e_reaxen = dft_energies['FS'] - dft_energies['IS']

    gap_co2form = gap_energies['CO2'] - gap_energies['CO'] - 0.5*gap_energies['O2']
    gap_oadfcc = gap_energies['O_FCC'] - gap_energies['Surf'] - 0.5*gap_energies['O2']
    gap_coadtop = gap_energies['CO_Top'] - gap_energies['Surf'] - gap_energies['CO']
    gap_fwbarr = gap_energies['TS_HCP'] - gap_energies['IS']
    gap_bkbarr = gap_energies['TS_HCP'] - gap_energies['FS']
    gap_reaxen = gap_energies['FS'] - gap_energies['IS']

    content = '===== Energy Statistics\n'
    content += ('{:<40s}{:>12.4f}{:>12.4f}\n').format('CO2 Formation Energy', e_co2form, gap_co2form)
    content += ('{:<40s}{:>12.4f}{:>12.4f}\n').format('O Adsorption Energy FCC', e_oadfcc, gap_oadfcc)
    content += ('{:<40s}{:>12.4f}{:>12.4f}\n').format('CO Adsorption Energy Top', e_coadtop, gap_coadtop)
    content += ('{:<40s}{:>12.4f}{:>12.4f}\n').format('Forward Reaction Barrier', e_fwbarr, gap_fwbarr)
    content += ('{:<40s}{:>12.4f}{:>12.4f}\n').format('Backward Reaction Barrier', e_bkbarr, gap_bkbarr)
    content += ('{:<40s}{:>12.4f}{:>12.4f}\n').format('Reaction Energy', e_reaxen, gap_reaxen)

    print(content)

    """
    # check lattice optimization
    atoms = train_frames[6].copy()
    atoms.set_calculator(gap)

    #print(atoms.get_stress(voigt=False))
    #print(atoms.get_calculator().get_virial(atoms))
    #print(atoms.get_potential_energy())
    print(atoms.get_cell())

    sf = StrainFilter(atoms)
    opt = BFGS(sf)
    opt.run(0.005)

    ase.io.write('opt.xyz',atoms)

    print(atoms.get_cell())

    # check surface formation energy
    atoms = train_frames[50].copy()
    atoms.set_calculator(gap)

    cons = FixAtoms(indices=[atom.index for atom in atoms if atom.position[2]<4.0])
    #print(atoms[0].position)
    atoms.set_constraint(cons)

    dyn = BFGS(atoms, trajectory='surf.traj')
    dyn.run(fmax=0.05)
    # -117.40321
    surf_traj = ase.io.read('./surf.traj', ':') # ase atom frames
    surf_opt = surf_traj[-1]
    ase.io.write('surf.xsd',surf_opt)
    """
