#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse

from ase.io import read, write
from ase.io.dftb import write_dftb

from vaspy.matstudio import XsdFile

def write_gen(atoms):
    poses = atoms.get_scaled_positions()
    
    return

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-x', '--xsd', \
            default='default.xsd', help='xsd File')
    parser.add_argument('-g', '--gen', \
            default='geo_in.gen', help='Gen Format File Name')

    args = parser.parse_args()

    # read xsd and get fixes
    if args.xsd == 'default.xsd':
        for fname in os.listdir('.'):
            name_list = fname.split('.')
            if name_list[-1] == 'xsd':
                xsdname = name_list[0]
                xsd = XsdFile(fname)
                break
    else:
        xsd = XsdFile(args.xsd)

    symbol_dict = {}
    for i, symbol in enumerate(xsd.atom_types):
        symbol_dict[symbol] = i+1

    symbols = []
    for num, sym in zip(xsd.atom_numbers,xsd.atom_types):
        symbols.extend([sym]*num)

    content = ('  {:<8d}  {:<s}\n').format(xsd.natom,'F')
    content += ('{:>4s}'*len(xsd.atom_types)+'\n').format(*xsd.atom_types)
    for i, (sym, coord) in enumerate(zip(symbols,xsd.data)):
        content += ('{:>8d}  {:>4d}  '+'{:>16.12f}  '*3+'\n')\
                .format(i+1,symbol_dict[sym],*coord)
    content += ('{:>16.12f}  '*3+'\n').format(*[0.]*3)
    for lat in xsd.bases:
        content += ('{:>16.12f}  '*3+'\n').format(*lat)

    with open(args.gen,'w') as writer:
        writer.write(content)

    fixes = xsd.tf
    moved_atoms = []
    for i, fix in enumerate(fixes):
        if fix[0] == 'T':
            moved_atoms.append(i+1)
    comment = 'MovedAtoms = {' + \
            ('{:<4d}  '*len(moved_atoms)).format(*moved_atoms) + \
            '}'
    print(comment)
