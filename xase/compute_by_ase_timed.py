#!/usr/bin/env python3

import os
import time
import json
import logging
import subprocess
import argparse

from pathlib import Path

import numpy as np

from ase.io import read, write
from ase.io.vasp import write_vasp 
from ase.calculators.vasp import Vasp2 # use Vasp in later ase version

# logger
logLevel = logging.INFO

logger = logging.getLogger(__name__)
logger.setLevel(logLevel)

fh = logging.FileHandler(filename='log.txt', mode='w')
fh.setLevel(logLevel)

ch = logging.StreamHandler()
ch.setLevel(logLevel)

logger.addHandler(ch)
logger.addHandler(fh)


basic_vasp_params = {
    # INCAR 
    "nwrite": 2, 
    "istart": 0, 
    "lcharg": False, 
    "lwave": False, 
    "lorbit": 10,
    "npar": 4,
    "xc": "pbe",
    "encut": 400,
    "prec": "Normal",
    "ediff": 1E-5,
    "nelm": 180, 
    "nelmin": 6, 
    "ispin": 1,
    "ismear": 1,
    "sigma": 0.2,
    "algo": "Fast", 
    "lreal": "Auto", 
    "isym": 0, 
    "ediffg": -0.05,
    "nsw": 200,
    "ibrion": 2,
    "isif": 2,
    "potim": 0.02, 
    # KPOINTS
    "gamma": True
}

VASP_COMMAND = os.environ['VASP_COMMAND']
VASP_PP_PATH = os.environ['VASP_PP_PATH']

def generate_calculator(calculator, kpts, command, directory='vasp-outputs'):
    """turn a dict into a ase-calculator"""
    calc_params = basic_vasp_params.copy()
    calc_params.update(kpts = kpts)
    calc_params.update(command = command)
    calc_params.update(directory = directory)

    calc = calculator(**calc_params)

    return calc

def create_vasp_inputs(directory, atoms, incar):
    # ===== initialise ===== 
    calc = Vasp2(
        command = VASP_COMMAND, 
        directory = 'dummy',  
        **basic_vasp_params, 
    )

    # convert to pathlib usage 
    if not directory.exists():
        directory.mkdir()
    else:
        pass

    calc.initialize(atoms)

    calc.string_params['system'] = Path(directory).name

    content = '\n>>>>> Modified ASE for VASP <<<<<\n'
    content += '    directory -> %s\n' %directory 
    #print(content)
    logger.info(content)

    # ===== POSCAR ===== 
    poscar = os.path.join(directory, 'POSCAR')
    write_vasp(poscar, atoms, direct=True, vasp5=True)

    # ===== KPOINTS =====
    kpts = np.linalg.norm(atoms.cell, axis=1).tolist()
    kpts = [int(20./k)+1 for k in kpts] 

    calc.input_params['gamma'] = True # use gamma-centred mesh
    calc.input_params['kpts'] = kpts 
    calc.write_kpoints(atoms=atoms, directory=directory)

    content = 'KPOINTS -->\n'
    content += "     Set k-point -> \033[1;35m{} {} {}\033[0m\n".format(*kpts)
    #print(content)
    logger.info(content)

    # ===== INCAR and 188 ===== 
    content = 'INCAR -->\n'

    calc.read_incar(incar)
    content += '    read template from %s\n' %incar 

    calc.write_incar(atoms, directory=directory)

    #print(content)
    logger.info(content)
    
    # ===== POTCAR =====
    calc.write_potcar(directory=directory)
    content = "POTCAR -->\n"
    content += "    Using POTCAR from %s\n" %VASP_PP_PATH
    #print(content)
    logger.info(content)

    # ===== VDW =====
    if calc.bool_params['luse_vdw']:
        calc.copy_vdw_kernel(directory=directory)

        content = "VDW -->\n"
        content += "    Using vdW Functional as %s\n" %calc.string_params['gga']
        #print(content)
        logger.info(content)

    return

def run_calculation(
        stru_file, indices, prefix, incar_list
    ):
    """"""
    # initialise few files
    for idx in range(len(incar_list)):
        with open('calculated_'+str(idx)+'.xyz', 'w') as writer:
            writer.write('')

    # read structures
    frames = read(stru_file, indices)
    start, end = indices.split(':')
    if start == '':
        start = 0
    else:
        start = int(start)
    logger.info('%d structures in %s from %d\n', len(frames), stru_file, start)

    # run calc
    for jdx, incar in enumerate(incar_list):
        logger.info(
            "\n===== Calculation Stage %d =====\n", jdx
        )
        for idx, atoms in enumerate(frames):
            logger.info(
                "\n===== Structure Number %d\n", idx+start
                #"Structure Number %d\n Index in XYZ is %d\n", idx, atoms.info['step']
            )
            
            directory = Path(prefix+'_'+str(jdx)+'_'+str(idx+start))
            create_vasp_inputs(directory, atoms, incar)

            try:
                proc = subprocess.run(
                    VASP_COMMAND, shell=True, cwd=directory, timeout=300,
                )

                if not proc.check_returncode():
                    vasprun = directory / 'vasprun.xml'
                    atoms = read(vasprun, format='vasp-xml')
                    atoms.info['source'] = directory.name
                    write('calculated_'+str(jdx)+'.xyz', atoms, append=True)
                    logger.info('successcode %s' %(proc.returncode))
                else:
                    logger.info('errorcode %s' %(proc.returncode))
            except subprocess.TimeoutExpired:
                logger.info('Time out...')
            except subprocess.CalledProcessError:
                logger.info('errorcode...' )

    return

if __name__ == '__main__':
    logger.info(
        '\nStart at %s\n', 
        time.asctime( time.localtime(time.time()) )
    )

    # args
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-s', '--structure', 
        default='example.xyz', 
        help='input structures stored in xyz format file'
    )
    parser.add_argument(
        '-p', '--parameter', 
        default='INCAR', nargs='*',
        help='calculator-related parameters in json format file'
    )
    parser.add_argument(
        '-i', '--indices', 
        default=':', 
        help='unsupported frame selection'
    )

    args = parser.parse_args()

    # run calculation 
    run_calculation(args.structure, args.indices, 'vasp', args.parameter)

    logger.info(
        '\nFinish at %s\n', 
        time.asctime( time.localtime(time.time()) )
    )
