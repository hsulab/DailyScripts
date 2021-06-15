#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import argparse

from ase.io import read, write

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-d', '--datadir', 
        help='time series files'
    )
    parser.add_argument(
        '-s', '--suffix', 
        default='data', help='time series files'
    )
    parser.add_argument(
        '-i', '--informat', 
        default='xsd', help='time series files'
    )
    parser.add_argument(
        '-o', '--outformat', 
        default='xsd', help='time series files'
    )

    args = parser.parse_args()

    stru_files = []
    frames = []
    data_dir = Path(args.datadir)
    for p in data_dir.glob('*.'+args.informat):
        stru_files.append(p)
    #stru_files = sorted(stru_files, key=lambda k: int(k.split('-')[-1]))
    stru_files.sort()

    for p in stru_files:
        print(p.stem + '.' + args.suffix)
        atoms = read(p)
        frames.append(atoms)
        write(data_dir / (p.stem + '.' + args.suffix), atoms, format=args.outformat)
    #write('merged.arc', frames, format='dmol-arc')
    write(data_dir / 'merged.xyz', frames)
