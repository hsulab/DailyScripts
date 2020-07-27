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
    gap = Potential(param_filename='./gap-2/GAP.xml')
    train_frames = ase.io.read('./train.xyz', ':') # ase atom frames
    print(len(train_frames))

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
    """
    surf_traj = ase.io.read('./surf.traj', ':') # ase atom frames
    surf_opt = surf_traj[-1]
    ase.io.write('surf.xsd',surf_opt)
