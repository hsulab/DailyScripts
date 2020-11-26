#!/usr/bin/env python3

import os
import time
import json
import logging
import argparse

from ase.io import read, write
from ase.calculators.vasp import Vasp2

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

# VASP environments

def read_structures(stru_file, indices):
    """"""
    frames = read(stru_file, ':')
    if type(frames) != list:
        frames = [frames]

    # TODO: add step info
    # for idx, atoms in enumerate(frames):
    #     atoms.info['step'] = idx

    logger.info('%d structures in %s\n', len(frames), stru_file)

    return frames

def parse_params(json_file):
    """pass json format params"""
    # read json
    with open(json_file, 'r') as fr:
        params = json.load(fr)

    # get calculator
    calc_prefix = 'calc'
    calc_command = params.pop('calc_command', None)
    if calc_command:
        if 'vasp' in calc_command:
            ase_calculator = Vasp2
            calc_options = params.pop('calc_options', None)
            if calc_options:
                pp_path = calc_options.pop('pp_path', None)
                if pp_path:
                    if 'VASP_PP_PATH' in os.environ.keys():
                        os.environ.pop('VASP_PP_PATH')
                    os.environ['VASP_PP_PATH'] = pp_path
                else:
                    raise ValueError(
                        'Not read key pp_path in calc_options in %s' %(json_file)
                    )
            else:
                raise ValueError(
                    'Not read key calc_options in %s' %(json_file)
                )
            logger.info('Use VASP to calculated structures.\n')
            calc_prefix = 'vasp'
        else:
            # TODO: add more calculators
            raise ValueError('Unsupported command in %s' %(json_file))
    else:
        raise ValueError('Not read key calc_command in %s' %(json_file))

    calc_params = params.pop('calc_params', None)
    if calc_params:
        if type(calc_params) == dict:
            calc_params = [calc_params]
        else:
            if type(calc_params) == list:
                pass
            else:
                raise ValueError('Cannot read calc_params in %s' %(json_file))
        logger.info('Total %d stages calculation.\n', len(calc_params))
    else:
        raise ValueError('Not read key calc_params in %s' %(json_file))

    return calc_prefix, ase_calculator, calc_command, calc_params

def generate_calculator(calculator, command, params, directory='vasp-outputs'):
    """turn a dict into a ase-calculator"""
    calc_params = params
    calc_params.update(command = command)
    calc_params.update(directory = directory)

    calc = calculator(**calc_params)

    return calc

def run_calculation(
        frames, prefix, calculator, calc_command, calc_params_list
    ):
    """"""
    # initialise few files
    for idx in range(len(calc_params_list)):
        with open('calculated_'+str(idx)+'.xyz', 'w') as writer:
            writer.write('')

    for jdx, calc_params in enumerate(calc_params_list):
        logger.info(
            "\n===== Calculation Stage %d =====\n", jdx
        )
        for idx, atoms in enumerate(frames):
            logger.info(
                "Structure Number %d\n", idx
                #"Structure Number %d\n Index in XYZ is %d\n", idx, atoms.info['step']
            )
            atoms.set_calculator(
                generate_calculator(
                    calculator, calc_command, calc_params, 
                    prefix+'_'+str(jdx)+'_'+str(idx)
                    #'vasp_'+str(jdx)+'_'+str(atoms.info['step'])
                )
            )
            dummy = atoms.get_forces() # just to call a calculation 
            write('calculated_'+str(jdx)+'.xyz', atoms, append=True)

        frames = read('calculated_'+str(jdx)+'.xyz', ':')
        if type(frames) != list:
            frames = [frames]

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
        default='params.json', 
        help='calculator-related parameters in json format file'
    )
    parser.add_argument(
        '-i', '--indices', 
        default=None, nargs='*', 
        help='unsupported frame selection'
    )

    args = parser.parse_args()

    # read structures
    frames = read_structures(args.structure, args.indices)

    # read calc params
    #calculator, calc_command, calc_params_list = parse_params('vasp_params.json')

    # run calculation 
    run_calculation(frames, *parse_params(args.parameter))

    logger.info(
        '\nFinish at %s\n', 
        time.asctime( time.localtime(time.time()) )
    )
