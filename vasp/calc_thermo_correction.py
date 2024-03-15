#!/usr/bin/env python

import os
import sys

import argparse
import subprocess
from pathlib import Path

from math import exp
from math import log

import numpy as np

def read_frequencies(outcar):
    command = "grep \"cm-1\" {}".format(outcar)
    proc = subprocess.Popen(
        command, shell=True, cwd=Path.cwd(),
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        encoding = "utf-8"
    )
    freqs = []
    lines = proc.stdout.readlines()
    for line in lines:
        data = line.strip().split()
        # print(data[1])
        if data[1] == "f":
            freq = complex(float(data[-2]), 0.)
        elif data[1] == "f/i=":
            freq = complex(0., float(data[-2]))
        else:
            raise ValueError("error in vasp freq...")
        freqs.append(freq)
    freqs = np.array(freqs)

    return freqs, lines

parser = argparse.ArgumentParser()
parser.add_argument(
    "outcar", default="./OUTCAR",
    help = "temperature for harmonic limit correction"
)
parser.add_argument(
    "-t", "--temperature", default=300., type=float,
    help = "temperature for harmonic limit correction"
)
parser.add_argument(
    "-ad", "--adsorbate", action="store_true",
    help = "whether calculate adsorbate"
)
args = parser.parse_args()

outcar = Path(args.outcar)
temperature = args.temperature
#read_frequencies("./VO-100/PGM2/IS/r1-freq/OUTCAR")
#freqs, lines = read_frequencies("./VO-100/PGM2/TS/r1-TS-freq/OUTCAR")
freqs, lines = read_frequencies(outcar)

use_ase = True
#use_ase = False
if use_ase:
    # get real vibrational energies
    vib_energies = [x for x in (np.real(freqs)/1000.).tolist() if x > 0.] # meV to eV

    # choose thermo
    is_adsorbate = args.adsorbate
    if is_adsorbate:
        from ase.thermochemistry import HarmonicThermo
        thermo = HarmonicThermo(vib_energies)
        print(thermo.get_helmholtz_energy(temperature))
    else:
        from ase.io import read
        contcar = outcar.parent / "CONTCAR"
        atoms = read(contcar)
        from ase.thermochemistry import IdealGasThermo
        thermo = IdealGasThermo(
            vib_energies=vib_energies,
            potentialenergy=0.,
            atoms=atoms,
            geometry="linear",
            symmetrynumber=2, 
            spin=0
        )
        print(thermo.get_gibbs_energy(temperature, pressure=101325.))
else:
    detail = True
    T = temperature
    strlist = lines
    wlist=[ float(str.split('2PiTHz')[1].split()[0]) for str in strlist if 'f/i' not in str ] # get cm^-1
    h = 6.63*10**(-34)
    c = 3*10**10
    KB = 1.38*10**(-23)
    Nmol = 6.023*10**23
    eV2J = 96485.310432
    R = 8.314
    hcw= [ h*c*x for x in wlist ]
    hcw_kbt= [ x/(KB*T) for x in hcw ]
    
    zp=sum([ x*Nmol/2 for x in hcw ]) #in J/mol
    S=sum([ R*(x/(exp(x)-1)-log(1-exp(-x))) for x in hcw_kbt ]) #in J/mol/K
    U=sum([ R*x/KB/(exp(x/KB/T)-1) for x in hcw ]) #in J/mol
    TS=T*S
    
    zp_ineV=zp/eV2J
    TS_ineV=TS/eV2J
    U_ineV=U/eV2J
    
    total = zp_ineV + U_ineV - TS_ineV
    print('Temperature = %s K\tZPE = %3.2f eV\tU = %3.2f eV\tTS = %3.2f eV\tZPE+U-TS = %3.2f eV\n'%(T,zp_ineV,U_ineV,TS_ineV,total))
    print(wlist)
