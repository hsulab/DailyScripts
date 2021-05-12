#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from ase.io import read, write 
from ase.calculators.emt import EMT 

if __name__ == '__main__':
    calc = EMT()
    evaluated_datfile = './evaluated.xyz'

    frames = read('./basic_strus.xyz', ':')
    for atoms in frames: 
        calc.reset() 
        atoms.calc = calc 
        dummy = atoms.get_forces() # single-point calculation 
        write(evaluated_datfile, atoms, append=True)

    pass
