#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import subprocess
import shutil

import copy

import numpy as np

from ase import units

from ase import Atoms
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
from ase.optimize import QuasiNewton
from ase.build import fcc111, add_adsorbate

# molecular dynamics
from ase.md.verlet import VelocityVerlet
from ase.md.nvtberendsen import NVTBerendsen
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution

import ase.io
import ase.io.dmol
from ase.io.trajectory import Trajectory


def write_plumed_dat(datPath, distance):
    content = "FLUSH STRIDE=1\n"
    content += "com1: COM ATOMS=1,2 MASS\n"
    content += "dis1: DISTANCE ATOMS=com1,16\n"
    content += "res1: RESTRAINT ARG=dis1 AT=%f KAPPA=200000.0\n" %distance
    content += "dis2: DISTANCE ATOMS=1,12\n"
    content += "dis3: DISTANCE ATOMS=1,8\n"
    content += "dis4: DISTANCE ATOMS=1,4\n"
    content += "lw1: LOWER_WALLS ARG=dis2 AT=0.30 KAPPA=100000.0\n"
    content += "lw2: LOWER_WALLS ARG=dis3 AT=0.30 KAPPA=100000.0\n"
    content += "lw3: LOWER_WALLS ARG=dis4 AT=0.30 KAPPA=100000.0\n"
    content += "PRINT STRIDE=1 ARG=dis1 FILE=colvar\n"
    content += "PRINT STRIDE=1 ARG=dis1,res1.bias FILE=bias\n"
    content += "PRINT STRIDE=1 ARG=dis2,lw1.bias,dis3,lw2.bias,dis4,lw3.bias FILE=walls\n"

    with open(datPath, 'w') as writer:
        writer.write(content)

    return

def generate_initial_structure(genPath, atoms_in, distance):
    distance = distance*10 # nm to angstrom
    # 'geo_in.gen'
    atoms = atoms_in.copy()
    
    CPos = atoms[0].position
    OPos = atoms[1].position
    PtPos = atoms[15].position

    CODistance = 1.12

    CPos = PtPos.copy()
    CPos[2] = PtPos[2] + distance - CODistance/2.

    OPos = PtPos.copy()
    OPos[2] = PtPos[2] + distance + CODistance/2.

    #print(atoms[1].position)
    #print(OPos)

    atoms.positions[0] = CPos.copy()
    atoms.positions[1] = OPos.copy()

    ase.io.write(genPath, atoms)

    return


if __name__ == '__main__':
    structureTemplatePath = './usTemp.xyz'
    atoms_in = ase.io.read(structureTemplatePath)

    #us_frames = ase.io.read('./us.xyz', ':')
    distances = np.linspace(0.29,0.33,5)
    #print(len(us_frames))
    for distance in distances:
        distance = round(distance, 2)
        directoryPath = 'd%4d' %(distance*10000)
        if os.path.exists(directoryPath):
            shutil.move(directoryPath, directoryPath+'.bak')
            os.makedirs(directoryPath)
            print(directoryPath)
        else:
            os.makedirs(directoryPath)
            print(directoryPath)
        shutil.copyfile('./dftb_in.hsd', directoryPath+'/dftb_in.hsd')
        genPath = os.path.join(directoryPath, 'geo_in.gen')
        generate_initial_structure(genPath, atoms_in, distance)
        #distance = np.linalg.norm(pos1 - pos2)
        write_plumed_dat(os.path.join(directoryPath, 'plumed.dat'), distance)
        #print(distance)
        #run_ase_plumed(atoms, potentialPath, directoryPath)
        #subprocess.Popen("echo 'haha' > dftb.out", shell=True, cwd=directoryPath)
        dftbPlus = "/home/mmm0586/apps/dftbplus/release/master/install/bin/dftb+"
        command = "mpirun -n 2 %s 2>&1 > dftb.out" %dftbPlus
        print(command)
        proc = subprocess.Popen(command, \
                shell=True, cwd=directoryPath)

        errorcode = proc.wait()
        if errorcode:
            path = os.path.abspath(self.directory)
            msg = ('Calculator "{}" failed with command "{}" failed in '
                   '{} with error code {}'.format('dftb+', command,
                                                  path, errorcode))
            print(msg)

