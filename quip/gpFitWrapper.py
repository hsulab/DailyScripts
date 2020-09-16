#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# gap_fit energy_parameter_name=energy force_parameter_name=forces do_copy_at_file=F sparse_separate_file=T gp_file=GAP.xml at_file=train.xyz default_sigma={0.008 0.04 0 0} gap={distance_2b cutoff=4.0 covariance_type=ard_se delta=0.5 theta_uniform=1.0 sparse_method=uniform add_species=T n_sparse=10}

# gap_fit \
#   at_file=train.xyz gp_file=GAP.xml \
#   do_copy_at_file=F sparse_separate_file=F \ 
#   energy_parameter_name=free_energy force_parameter_name=forces virial_parameter_name=virial \ 
#   default_sigma={0.004 0.04 0.04 0.0} \
#   gap={ \
#       distance_2b cutoff=4.0 delta=0.8 theta_uniform=1.0 \ 
#           covariance_type=ard_se sparse_method=uniform n_sparse=100 : \
#       angle_3b cutoff=4.0 delta=0.4 theta_uniform=1.0 \ 
#           covariance_type=ard_se sparse_method=uniform n_sparse=250 : \
#       soap cutoff=4.0 delta=0.2 zeta=2 atom_sigma=0.7 l_max=6 n_max=6 \ 
#           covariance_type=dot_product sparse_method=cur_points n_sparse=250 \
#       }

import os
import subprocess

import argparse

import ase.io

"""
This is a gap_fit wrapper for quickly generating the gap_fit command.
"""

file_keys = [
        'at_file', # training data file
        'gp_file', # output gp file
        ]

bool_keys = [
        'do_copy_at_file', # 
        'sparse_separate_file', # sparse points file
        ]

string_keys = [
    'energy_parameter_name',
    'force_parameter_name',
    'virial_parameter_name',
    'hessian_parameter_name',
    ]

list_keys = [
        'default_sigma'
        ]

descriptor_types = [
        'distance_2b', 
        'angle_3b',
        'soap'
        ]


class GapFit(object):
    """
    """

    #executable = 'gerun gap_fit_mpi'
    executable = 'gap_fit_mpi'

    def __init__(self, directory='.', executable=executable):
        # energy force virial hessian
        self.energy_parameter_name = 'energy'

        # gap_fit
        self.executable = executable

        # 
        self.directory = directory
        #if os.path.exists(self.directory):
        #    raise ValueError('Please check directory %s.' %self.directory)

        # files related options
        self.at_file = 'train.xyz'
        self.gp_file = 'GAP.xml'

        self.do_copy_at_file = 'F'
        self.sparse_separate_file = 'F'

        # default_sigma
        self.default_sigma = [0.,0.,0.,0.]

        # 
        self.para_names = {}
        self.gaps = []

    def add_para_names(self, **kwargs):
        para_args = [
            'energy_parameter_name',
            'force_parameter_name',
            'virial_parameter_name',
            'hessian_parameter_name',
            ]
        para_keys = [p.split('_')[0] for p in para_args]

        for key, value in kwargs.items():
            if key in para_keys:
                idx = para_keys.index(key)
                self.para_names.update({para_args[idx]: value})
            else:
                raise ValueError('Key must be in'+('{}'*len(para_keys)).format(*para_keys))

        return

    def set_default_sigma():
        pass

    def add_gap_command(self, desc_list):
        """
        """
        self.gap_command = ':'.join(desc_list)

        return 

    def generate_command(self):
        # first add exe
        command = self.executable + ' \\\n'

        # add parameter_name
        if self.para_names:
            para_command = []
            for key, value in self.para_names.items():
                para_command.append(key + '=' + value)
            para_command = ' '.join(para_command)
        else:
            raise ValueError('No parameter_names')
        command += para_command + ' \\\n'

        # add balabla
        command += ' do_copy_at_file=F sparse_separate_file=F gp_file=GAP.xml at_file=train.xyz default_sigma={0.008 0.04 0.0 0} '

        command += 'gap={'
        command += self.gap_command
        command += '}'

        return command

    def run(self):
        command = self.generate_command()

        try:
            proc = subprocess.Popen(command, shell=True, cwd=self.directory)
        except OSError as err:
            msg = 'Failed to execute "{}"'.format(command)
            raise EnvironmentError(msg) from err

        errorcode = proc.wait() # wait for child process to terminate

        if errorcode:
            path = os.path.abspath(self.directory)
            msg = ('Calculator "{}" failed with command "{}" failed in '
                '{} with error code {}'.format(self.name, command,
                                                path, errorcode))
            raise ValueError(msg)

if __name__ == '__main__':
    # parser
    parser = argparse.ArgumentParser()

    parser.add_argument('-m', '--mode', \
            default='com', help='mode')
    parser.add_argument('-d', '--dat', \
            default='sampled_structures.xyz', help='datafile in total')

    args = parser.parse_args()

    # descriptors
    import gpFitDescriptors as gpd

    #
    gf = GapFit(directory='../training/gap-5')

    gf.add_gap_command(gpd.desc_list)

    # para names
    gf.add_para_names(energy='free_energy',force='forces',virial='haha')

    # k-fold
    total_frames = ase.io.read(args.dat, ':')
    from sklearn.model_selection import KFold
    kf = KFold(n_splits=5, shuffle=True, random_state=42)
    for i in kf.split(total_frames):
        indices = i
        break
    train_indices = indices[0]
    valid_indices = indices[1]

    train_frames, valid_frames = [], []
    for idx in train_indices:
        train_frames.append(total_frames[idx])
    for idx in valid_indices:
        valid_frames.append(total_frames[idx])

    ase.io.write('train.xyz', train_frames)
    ase.io.write('valid.xyz', valid_frames)

    if args.mode == 'run':
        gf.run()
    else:
        print(gf.generate_command())

