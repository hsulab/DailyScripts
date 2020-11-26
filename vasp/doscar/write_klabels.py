#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from numpy import pi
import numpy as np 


def real2recip(poscar='POSCAR'):
    # read POSCAR
    with open(poscar, 'r') as reader:
        lines = reader.readlines()
        lines = [line.strip('\n').split() for line in lines]

        scale = float(lines[1][0])
        real_lat = np.array(lines[2:5], 'float')*scale # real lattice vector

        a1, a2, a3 = real_lat[0], real_lat[1], real_lat[2]
        volume = np.dot(a1, np.cross(a2, a3))

        b1 = 2*pi*np.cross(a2, a3) / volume
        b2 = 2*pi*np.cross(a3, a1) / volume
        b3 = 2*pi*np.cross(a1, a2) / volume

        recip_lat = np.array([b1, b2, b3])

        return recip_lat


def write_klabels():
    with open('KPATH.in', 'r') as reader:
        lines = reader.readlines()
    lines = [line.strip('\n').split() for line in lines]

    kpaths = []
    for line in lines[4:]:
        if len(line) == 4:
            kpaths.append([line[3], np.array(line[:3], 'float')])

    recip_lat = real2recip()
    klabels, kcoords = [], []
    for i, (label, coord) in enumerate(kpaths):
        if i%2 == 0:
            if i == 0:
                klabels.append(label)
                kcoords.append(coord)
            else:
                last_label = kpaths[i-1][0]
                if last_label == label:
                    klabels.append(label)
                    kcoords.append(coord)
                else:
                    klabels.append(last_label + '|' + label)
                    kcoords.append(kpaths[i-1][1])
        elif i%2 == 1 and i+1 == len(kpaths):
            klabels.append(label)
            kcoords.append(coord)
    
    ksteps = []
    content = '%10s%50s\n' %('K-Label', 'K-Coordinate')
    for i, (label, coord) in enumerate(zip(klabels, kcoords)):
        if i == 0:
            kstep = np.linalg.norm(np.dot(coord, recip_lat))
            content += '%10s%15.3f\n' %(label, kstep)
            ksteps.append(kstep)
        else:
            last_coord = kcoords[i-1]
            kstep = ksteps[i-1] + \
                    np.linalg.norm(np.dot(coord-last_coord, recip_lat))
            content += '%10s%15.3f\n' %(label, kstep)
            ksteps.append(kstep)
    print(content)


if __name__ == '__main__':
    write_klabels()
