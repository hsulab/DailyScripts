#!/usr/bin/env python
'''
Notes:
    Script to convert .xsd file to VASP input files.
    This version is specialized for Thomas the Server.
    Last updated by jyxu.
Description:
    Create VASP input files, POSCAR->POTCAR->PBS.sh->INCAR->...
'''

import argparse
import os
import re
import sys
import logging

import numpy as np

from vaspy.matstudio import XsdFile
from vaspy.incar import InCar, Param
from vaspy.vasp_para_db import *
from vaspy import PY2

if PY2:
    import commands as subprocess
else:
    import subprocess

_logger = logging.getLogger("vaspy.script")

POTDIR = r'/home/mmm0586/apps/vasp/pot/potpaw_PBE/'
PBSSCRIPT = r'vasp.script'


def read_xsd():
    # Get xsd file
    for filename in os.listdir('.'):
        name_list = filename.split('.')
        if name_list[-1] == 'xsd':
            xsdname = name_list[0]
            xsd = XsdFile(filename)
            break

    content = 'XSD -->\n'
    content += "    The xsdfile is : \033[31;1;4m%s\033[0m \n" %filename

    return content, xsd, xsdname

def create_poscar(xsd):
    """"""
    # read xsd and get poscar content
    poscar_content = xsd.get_poscar_content(bases_const=1.0)
    with open('POSCAR', 'w') as f:
        f.write(poscar_content)
    
    # get xsd info
    elements, numbers = xsd.atom_types, xsd.atom_numbers
    symbols = []
    for e, n in zip(elements, numbers):
        symbols.extend([e]*n)

    tf_dict = {}
    for e in elements:
        tf_dict[e] = [0, 0]

    for s, tf in zip(symbols, xsd.tf):
        if 'T' in tf:
            tf_dict[s][0] += 1
        elif 'F' in tf:
            tf_dict[s][1] += 1

    content = '\nPOSCAR -->\n'
    content += '    \033[4m Element       Numbers      \033[0m\n'
    t_count, f_count = 0, 0
    for s, tf in tf_dict.items():
        content += '         %2s \033[1;33m%4d\033[0m \033[32m%4d\033[0m(T)\033[31m%4d\033[0m(F)\n'\
                %(s, tf[0]+tf[1], tf[0], tf[1])
        t_count += tf[0]
        f_count += tf[1]
    content += '    \033[4m                            \033[0m\n'
    content += '      %2s \033[1;33m%4d\033[0m \033[32m%4d\033[0m(T)\033[31m%4d\033[0m(F)\n'\
                %('Total', t_count+f_count, t_count, f_count)
    content += '\n'

    return content

def create_potcar(xsd):
    potdir = POTDIR
    # delete old POTCAR
    if os.path.exists('./POTCAR'):
        os.remove('./POTCAR')

    for elem in xsd.atom_types:
        if os.path.exists(potdir + elem):
            potcar = potdir + elem + '/POTCAR_' + elem
        else:
            print('No POTCAR for ' + elem)
            sys.exit(1)
        subprocess.getstatusoutput('cat ' + potcar + ' >> ./POTCAR')

    content = 'POTCAR -->\n'
    content += ('     POTCAR Merged with \033[1;35m '\
            + '{:<3s}'*len(xsd.atom_types) + ' \033[0m\n').format(*xsd.atom_types)
    content += '\n'

    return content

def create_kpoints(args, xsd):
    # create kpoints
    if not args.kpoints:
        kpoints = []
        for base in xsd.bases:
            l = np.dot(base, base)**0.5
            kpt = int(20/l) + 1
            kpoints.append(str(kpt))
    else:
        kpoints = [i.strip() for i in args.kpoints.split(",")]

    kpt_str = ' '.join(kpoints)
    kpt_content = 'mesh auto\n0\nG\n' + kpt_str + '\n0 0 0\n'

    with open('KPOINTS', 'w') as f:
        f.write(kpt_content)

    content = 'KPOINTS -->\n'
    content += "     Set k-point -> \033[1;35m{} {} {}\033[0m\n".format(*kpoints)
    content += '\n'

    return content

def create_pbs(xsdname, args):
    jobname = xsdname
    with open('vasp.script', 'r') as f:
        content_list = f.readlines()

    content = 'PBS -->\n'

    pbs_content = ''
    for i, line in enumerate(content_list):
        if line.startswith('#$ -N'):
            pbs_content += "#$ -N %s\n" %jobname
            content += "     job name -> \033[1;35m{}\033[m\n".format(jobname)
        elif line.startswith('#$ -pe'):
            if args.ncpu:
                ncpu = args.ncpu
                pbs_content += '#$ -pe mpi %s\n' %args.ncpu
            else:
                ncpu = line.strip().split()[-1]
                pbs_content += line
            content += "     Number of Cores -> \033[1;35m{}\033[0m ".format(ncpu)
        elif line.startswith('#$ -P'):
            if args.queue:
                queue = args.queue
                pbs_content += '#$ -P %s\n' %args.queue
            else:
                queue = line.strip().split()[-1]
                pbs_content += line
            content += "     Queue -> \033[1;35m{}\033[0m\n".format(queue)
        else:
            pbs_content += line

    content += '\n'

    with open(PBSSCRIPT, 'w') as f:
        f.writelines(pbs_content)

    return content

def create_incar(xsdname, args):
    """"""
    if os.path.exists('INCAR'):
        return 'INCAR already exists, not changed.\n'

    incar = InCar()
    if args.task:
        task = args.task
        incar.quickgen(args.task)
    else:
        task = 'SC'
    content = 'INCAR -->\n'
    content += '     Task is %s\n' %task

    if PY2:
        pname_value_pairs = args.__dict__.iteritems()
    else:
        pname_value_pairs = args.__dict__.items()

    for pname, value in pname_value_pairs :
        if (value is not None):
            if pname in incar.pnames:
                incar.set(pname, value)
            else:
                if (pname in INCAR_PARAMETERS.keys()):
                    incar.add(pname, value)
            content += "     {} -> {}\n".format(pname, value)
    #incar.set('SYSTEM', xsdname)

    incar.tofile(verbos=False)

    content += '\n'

    return content

def create_ts(xsd):
    """"""
    atom_idxs, atom_names = [], []
    for idx, atom_name in enumerate(xsd.atom_names):
        if atom_name.endswith('_c'):
            atom_idxs.append(idx)
            atom_names.append(atom_name)
    # If constrained get distance and create fort.188
    if atom_idxs:
        if len(atom_idxs) > 2:
            raise ValueError("More than two atoms end with '_c'")
        pt1, pt2 = [xsd.data[idx, :] for idx in atom_idxs]

        # Use Ax=b convert to cartisan coordinate
        diff = pt1 - pt2
        A = np.matrix(xsd.bases.T)
        x = np.matrix(diff).T
        b = A*x
        distance = np.linalg.norm(b)

        # Create fort.188
        ts_content = '1\n3\n6\n4\n0.04\n%-5d%-5d%f\n0\n' % \
            (atom_idxs[0]+1, atom_idxs[1]+1, distance)

        with open('fort.188', 'w') as f:
            f.write(ts_content)

        content = 'Constrain TS -->\n'
        content += "     fort.188 has been created.\n"
        content += '     ' + '-'*20 + '\n'
        content += "     atom number: {:<5d}{:<5d}\n".format(atom_idxs[0]+1, atom_idxs[1]+1)
        content += "     atom name: {} {}\n".format(*atom_names)
        content += "     distance: {:f}\n".format(distance)
        content += '     ' + '-'*20 + '\n'

        # Set IBRION = 1
        #incar.set('IBRION', 1)
        #content += "{} -> {}\n".format("IBRION", "1")
        content += '     Note: IBRION must be 1.\n'

        content += '\n'

        return content
    else:
        return ''


if "__main__" == __name__:
    # Copy vasp.script
    # subprocess.getstatusoutput('cp $HOME/example/vasp.script ./')

    # Set argument parser.
    parser = argparse.ArgumentParser()

    # Add optional argument.
    parser.add_argument("-t", "--task", nargs='?', const='SC', help="set task type")
    parser.add_argument("-k", "--kpoints", help="set k-points")
    parser.add_argument("-q", "--queue", choices=('Gold', 'Test'), help="pbs queue type")
    parser.add_argument("-n", "--ncpu", help="cpu number in total")

    # Add all possible arguments in INCAR file.
    parameters = INCAR_PARAMETERS.keys()
    for parameter in parameters:
        help_info = "%s" %INCAR_PARAMETERS[parameter][2]
        parser.add_argument("--{}".format(parameter), help=help_info)

    # add all parses
    args = parser.parse_args()
    print('{:<20}{:^20}{:<20}'\
            .format('>'*20, 'VASPy Info', '<'*20))

    # start create vasp inputs
    content = '{:<20}{:^20}{:<20}\n'\
            .format('>'*20, 'Check VASP Inputs', '<'*20)

    # read xsdfile
    c, xsd, xsdname = read_xsd()
    content += c
    
    # Create POSCAR
    c = create_poscar(xsd)
    content += c

    # Create POTCAR
    c = create_potcar(xsd)
    content += c
    
    # Creat KPOINTS
    c = create_kpoints(args, xsd)
    content += c

    # Get content line list.
    #c = create_pbs(xsdname, args)
    #content += c
    
    # Create fort.188
    #c = create_ts(xsd)
    #content += c

    # set and write INCAR
    #c = create_incar(xsdname, args)
    #content += c
    
    content += '{:<20}{:^20}{:<20}\n'\
            .format('>'*20, 'End Check', '<'*20)
    print('{:<20}{:^20}{:<20}'\
            .format('>'*20, 'END VASPy', '<'*20))
    
    print(content)
