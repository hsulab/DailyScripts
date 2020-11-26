#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from numpy import pi
import numpy as np

"""
Author: jyxu
Description:
    A reproduction of band structure process of vaspkit in Python, 
    which avoids a malloc bug due to a large amount of bands.
Notes:
    Not compatible with spin calculation well.
"""


def real2recip(poscar='POSCAR'):
    """Read POSCAR and calculate reciprocal lattice."""
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
    #print(recip_lat)

    return recip_lat


def read_fermi():
    """Read DOSCAR and get fermi-energy."""
    # read doscar
    with open('DOSCAR', 'r') as reader:
        lines = reader.readlines()
    lines = [line.strip('\n').split() for line in lines]

    # read 
    e_max, e_min, nrows, e_fermi = [float(i) for i in lines[5][:4]]

    return e_max, e_min, e_fermi


def read_eigenvalue_info():
    """Only read first 6 lines of EIGENVAL."""
    # read eigenvalue
    with open('EIGENVAL', 'r') as reader:
        lines = reader.readlines()
    lines = [line.strip('\n').split() for line in lines]

    spin = int(lines[0][-1])
    nelectrons, nkpt, nband = [int(i) for i in lines[5]]

    return lines, nelectrons, nkpt, nband, spin


def read_eigenvalue(lines, nelectrons, nkpt, nband, spin):
    """
    Read EIGENVAL and get band energy.
    if spin==1, read first col elseif spin==2, read second col.
    """
    # read eigenvalue

    # get level, kmesh
    kmesh, level = [], []
    for i in range(nkpt):
        start_index = i*(nband + 2) + 6
        kmesh.append(lines[start_index+1][:3])

        if spin == 1 or spin == 2:
            klevel = [col[spin] for col in lines[start_index+2:start_index+2+nband]]
            level.append(klevel)

    kmesh = np.array(kmesh, 'float64') # kpoints
    level = np.array(level, 'float64') # band energies

    # get reciprocal lattice
    recip_lattice = real2recip()

    # get kstep
    kstep = np.zeros((nkpt,1))
    for i in range(nkpt-1):
        kspc = kmesh[i+1] - kmesh[i] # vector in reciprocal space
        kspc = np.dot(kspc, recip_lattice) # vector in real space
        kstep[i+1] = kstep[i] + np.linalg.norm(kspc)

    return kmesh, level, kstep


def read_gap():
    """Calculate gap using SPIN-UP bands."""
    lines, nelectrons, nkpt, nband, spin = read_eigenvalue_info()
    kmesh, level, kstep = read_eigenvalue(lines, nelectrons, nkpt, nband, 1)
    
    if nelectrons%2 == 0:
        nth_band_vbm = int(nelectrons / 2)
        nth_band_cbm = int(nelectrons / 2 + 1)
    else:
        raise ValueError('Band for Odd Electrons is Not Supported.')
    #print(nth_band_vbm, nth_band_cbm)

    # pay attetion to the index!!!
    vbm = np.max(level[:,nth_band_vbm-1])
    #print(level[:,nth_band_vbm])
    cbm = np.min(level[:,nth_band_cbm-1])
    #print(level[:,nth_band_cbm])

    for ikpt in range(nkpt):
        if level[ikpt, nth_band_vbm-1] >= vbm:
            vbm = level[ikpt, nth_band_vbm-1]
            vbm_kpt_coord = kmesh[ikpt]
        if level[ikpt, nth_band_cbm-1] <= cbm:
            cbm = level[ikpt, nth_band_cbm-1]
            cbm_kpt_coord = kmesh[ikpt]
    #print(vbm_kpt_coord, cbm_kpt_coord)

    band_gap = cbm - vbm 

    # calc band character
    if band_gap < 1e-6:
        band_character = 'Metallic'
    else:
        if np.linalg.norm(vbm_kpt_coord - cbm_kpt_coord) > 1e-6:
            band_character = 'Indirect'
        else:
            band_character = 'Direct'

    return band_character, band_gap, vbm, cbm, \
            nth_band_vbm, nth_band_cbm, vbm_kpt_coord, cbm_kpt_coord,


def write_gap():
    # read gap
    band_character, band_gap, vbm, cbm, \
            nth_band_vbm, nth_band_cbm, vbm_kpt_coord, cbm_kpt_coord = read_gap()

    content = '--BAND GAP--\n'
    content += '%25s%10s\n' %('Band Character:', band_character)
    content += '%25s%10.4f\n' %('Band Gap (eV):', band_gap)
    content += '%25s%10.4f\n' %('Eigenvalue of VBM (eV):', vbm)
    content += '%25s%10.4f\n' %('Eigenvalue of CBM (eV):', cbm)
    content += '%25s%10.4f\n' %('Fermi Energy (eV):', 0)
    content += '%25s%10d%10d\n' \
            %('HOMO and LUMO bands:', nth_band_vbm, nth_band_cbm)
    content += '%25s%10.6f%10.6f%10.6f\n' %('Location of VBM:', \
            vbm_kpt_coord[0], vbm_kpt_coord[1], vbm_kpt_coord[2])
    content += '%25s%10.6f%10.6f%10.6f\n' %('Location of CBM:', \
            cbm_kpt_coord[0], cbm_kpt_coord[1], cbm_kpt_coord[2])
    content += '--END--\n'
    
    # write gap 
    with open('BAND_GAP', 'w') as writer:
        writer.write(content)


def write_single_band(band_dat, nband, kstep, level, e_fermi):
    # read KPOINTS
    with open('KPOINTS', 'r') as reader:
        lines = reader.readlines()
    line = lines[0].strip('\n').split(':')[1].split()
    nkpt_ibz, ndiv, bandpath_num = [int(i) for i in line[2:5]]

    nkpt = ndiv # number of kpoints with zero-weight

    kstep_0 = kstep[nkpt_ibz]
    kstep = kstep[nkpt_ibz:] - kstep_0

    level = level[nkpt_ibz:] - e_fermi
    
    '''
    for iband in range(nband):
        for ikpt in range(nkpt):
            level[ikpt,iband] -= e_fermi 
            if level[ikpt,iband] > e_max:
                level[ikpt,iband] = e_max
            if level[ikpt,iband] < e_min:
                level[ikpt,iband] = e_min
    '''
    
    # write reformatted.dat
    content = '#KPATH\n'
    for ikpt in range(nkpt):
        content += '%6.2f' %kstep[ikpt]
        for iband in range(nband):
            content += '%9.3f' %level[ikpt][iband]
        content += '\n'

    with open(band_dat, 'w') as writer:
        writer.write(content)


def write_band():
    # read DOSCAR
    e_max, e_min, e_fermi = read_fermi()

    # read band
    lines, nelectrons, nkpt, nband, spin = read_eigenvalue_info()
    if spin == 1:
        band_dat = 'REFORMATTED_BAND.dat'
        kmesh, level, kstep = read_eigenvalue(lines, nelectrons, nkpt, nband, 1)

        # write band.dat
        content ='#KPATH\n'
        for iband in range(nband):
            content += '#Band-index%5d\n' %(iband+1)
            for ikpt in range(nkpt):
                if iband%2 == 0:
                    content += '%14.5f%14.5f\n' \
                            %(kstep[ikpt], level[ikpt][iband])
                elif iband%2 == 1:
                    content += '%14.5f%14.5f\n' \
                            %(kstep[nkpt-ikpt-1], level[nkpt-ikpt-1][iband])

        with open('BAND.dat', 'w') as writer:
            writer.write(content)

        write_single_band(band_dat, nband, kstep, level, e_fermi)
    elif spin == 2:
        band_up, band_dw = 'REFORMATTED_BAND_UP.dat', 'REFORMATTED_BAND_DW.dat'
        kmesh, level, kstep = read_eigenvalue(lines, nelectrons, nkpt, nband, 1)
        write_single_band(band_up, nband, kstep, level, e_fermi)
        kmesh, level, kstep = read_eigenvalue(lines, nelectrons, nkpt, nband, 2)
        write_single_band(band_dw, nband, kstep, level, e_fermi)


def write_klabels():
    with open('KPATH.in', 'r') as reader:
        lines = reader.readlines()
    lines = [line.strip('\n').split() for line in lines]
    
    # high symmetry k-points
    symkpoints, kpoints = {}, []
    for line in lines[4:]:
        if len(line) == 4:
            symkpoints[line[3]] = np.array(line[:3], 'float')
            kpoints.append([line[3], np.array(line[:3], 'float')])

    # get kpaths
    kpath, kpaths = [], [] 
    for i, (label, coord) in enumerate(kpoints):
        if i == 0:
            kpath.append(label)
        else:
            last_label = kpoints[i-1][0]
            if i%2 == 0:
                if label == last_label:
                    kpath.append(label)
                else:
                    kpath_new = kpath.copy()
                    kpath_new.append(last_label)
                    kpaths.append(kpath_new)
                    kpath.append(label)
            elif i%2 == 1:
                if i+1 == len(kpoints):
                    kpath.append(label)
                    kpaths.append(kpath)

    # add coordinate in path
    recip_lat = real2recip()
    ksteps = []
    for kpath in kpaths:
        kstep = []
        for i, label in enumerate(kpath):
            if i == 0:
                kstep.append(np.linalg.norm(symkpoints[label]))
            else:
                last_label = kpath[i-1]
                kspc = np.linalg.norm(np.dot(\
                        symkpoints[label] - symkpoints[last_label], recip_lat))
                kstep.append(kstep[i-1] + kspc)
            #print('%10s%15.3f' %(label, kstep[i]))
        #print('')
        ksteps.append(kstep)

    # write KLABELS
    content = 'K-Label    K-Coordinate in band-structure plots\n'
    for i in range(len(kpaths)):
        kpath, kstep = kpaths[i], ksteps[i]
        if i == 0: 
            for j in range(len(kpath)):
                label, coord = kpath[j], kstep[j]
                if len(kpaths) != i+1 and j+1 == len(kpath):
                    label = label + '|' + kpaths[i+1][j]
                #print('%10s%15.3f' %(label, coord))
                content += '%10s%15.3f\n' %(label, coord)
        elif i+1:# != len(kpaths):
            coord = last_coord
            for j in range(len(kpaths[i-1]), len(kpath)):
                label = kpath[j]
                coord += kstep[j] - kstep[j-1]
                if len(kpaths) != i+1 and j+1 == len(kpath):
                    label = label + '|' + kpaths[i+1][j]
                content += '%10s%15.3f\n' %(label, coord)
        last_coord = coord

    with open('KLABELS', 'w') as writer:
        writer.write(content)


if __name__ == '__main__':
    #read_fermi()
    #read_eigenvalue()
    #real2recip()
    #read_gap()
    write_gap()
    write_band()
    write_klabels()
