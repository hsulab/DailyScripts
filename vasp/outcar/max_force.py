#!/usr/bin/env python

import argparse
import logging
import os
import sys

from vaspy import PY2
if PY2:
    import commands as subprocess
else:
    import subprocess

import numpy as np

from vaspy.incar import InCar
from vaspy.atomco import PosCar
from vaspy.iter import OutCar
from vaspy.matstudio import XsdFile
from vaspy.functions import str2list

_logger = logging.getLogger("vaspy.script")

# Set arguments parser.
parser = argparse.ArgumentParser()
parser.add_argument("--xsd", help="create MaterStudio .xsd file",
                    action="store_true")
args = parser.parse_args()

outcar = OutCar()
poses, forces = outcar.forces(-1)

poscar = PosCar()
tfs = poscar.tf

incar = InCar()
#print(incar.paras)
EDIFFG = incar.get_pvalue('EDIFFG')
#EDIFFG = 0.05

bad_ids = []
max_idx = [-1, -1, -1]
max_force = [0., 0., 0.]
for i, (pos, force, tf) in enumerate(zip(poses, forces, tfs)):
    for j, (c, f) in enumerate(zip(tf, force)):
        if c == 'T':
            if np.fabs(f) >= max_force[j]:
                max_idx[j] = i
                max_force[j] = np.fabs(f)
            if np.fabs(f) > np.fabs(float(EDIFFG)):
                if i not in bad_ids:
                    bad_ids.append(i)

for i, x in zip(max_idx, ['x', 'y', 'z']):
    if i != -1:
        _logger.info("{:<15s}: {}".format("max force atom in %s" %x, i+1))
        _logger.info("{:<15s}: ({}, {}, {})".format("atom position", *poses[i]))
        _logger.info("{:<15s}: {}, {}, {}".format("forces", *forces[i]))
    else:
        _logger.info("Cannot get the max force in %s." %x)

_logger.info("Atoms not optimized >>> \n"+" "*24+"%s" %([i+1 for i in bad_ids]))

# Get fort.188 info.
if os.path.exists('./fort.188'):
    with open('fort.188', 'r') as f:
        atom_info = f.readlines()[5]
    _logger.info("{:<10s}{:<10s}{:<15s}".format("Atom1", "Atom2", "DISTANCE"))
    _logger.info("-"*30)
    _logger.info("{:<10s}{:<10s}{:<15s}\n".format(*str2list(atom_info)))

# Create .xsd file.
if args.xsd:
    status, output = subprocess.getstatusoutput('ls *.xsd | head -1')
    if not output.endswith('.xsd'):
        _logger.info("No .xsd file in current directory.")
        sys.exit(1)
    xsd = XsdFile(filename=output)
    # modify atom color
    xsd.modify_color(atom_number=outcar.last_max_atom)
    jobname = output.split('.')[0]
    filename = jobname + '-force.xsd'
    xsd.tofile(filename=filename)
    _logger.info(filename + " has been created.")

