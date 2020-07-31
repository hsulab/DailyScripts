#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os

import numpy as np

import ase.io
from ase import Atoms, Atom


def write_xyz(symbols,positions,**kwargs):
    """positions in cartesian (AA) and forces in eV/AA"""
    # check number of atoms
    natoms = len(symbols)

    # implented properties
    comment_properties = [\
            'step',\
            'pbc','Lattice',\
            'energy','free_energy',\
            'stress','virial']
    atomic_properties = ['species','pos','forces']

    # check args
    comment_content = ''
    if len(kwargs) > 0:
        for key in comment_properties:
            if key in kwargs.keys():
                value = kwargs[key]
                if key in ['energy','free_energy']:
                    # float
                    value = float(value)
                    comment_content += ("{:<s}="+"{:<.4f}"+" ") \
                            .format(key,value)
                elif key in ['Lattice','stress','virial']:
                    # list of float properties
                    value = list(np.array(value,dtype=float).ravel())
                    if len(value) != 9:
                        raise ValueError('Lattice/stress/virial must have 9 components.')
                    comment_content += ("{:<s}="+"\""+"{:<.4f} "*len(value)+"\""+" ") \
                            .format(key,*value)
                elif key in ['pbc']:
                    # list of strings
                    comment_content += ("{:<s}="+"\""+"{:<s} "*len(value)+"\""+" ") \
                            .format(key,*value)
                elif key in ['step']:
                    value = int(value)
                    comment_content += ("{:<s}="+"{:<d} "+" ") \
                            .format(key,value)
            #else:
            #    raise ValueError('Unsupported properties in extended-xyz.')
    else:
        pass

    # write content
    content = "{:<d}\n".format(natoms)

    if 'forces' in kwargs.keys():
        forces = kwargs['forces']
        comment_content += 'Properties=species:S:1:pos:R:3:forces:R:3\n'
        content += comment_content
        for i in range(natoms):
            content += ('{:<4s} '+'{:>12.6f} '*6+'\n')\
                .format(symbols[i],*list(positions[i]),*list(forces[i]))
    else:
        comment_content += 'Properties=species:S:1:pos:R:3\n'
        content += comment_content
        for i in range(natoms):
            content += ('{:<4s} '+'{:>12.6f} '*3+'\n')\
                .format(symbols[i],*list(positions[i]))

    return content


if __name__ == '__main__':
    pass
