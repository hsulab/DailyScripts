#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse

def read_charges(datfile):
    with open(datfile, 'r') as reader:
        lines = reader.readlines()
    lines = [line.strip().split() for line in lines]
    lines = lines[2:-4]

    charges = []
    for line in lines:
        charges.append(float(line[4]))

    return charges

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-f', '--datafile')
    parser.add_argument('-f2', '--datafile2')

    args = parser.parse_args()

    charges1 = read_charges(args.datafile)
    charges2 = read_charges(args.datafile2)

    pt = 0
    for i, (c1, c2) in enumerate(zip(charges1, charges2)):
        if i >= 72 and i < 83:
            #print(c1-c2)
            pt += c1-c2

    print(pt)
