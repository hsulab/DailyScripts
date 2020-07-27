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

    executable = 'gap_fit'

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

    def add_gap(self, desc, spar):
        """
        desc = {
        'name': 'distance_2b',
        'cutoff': 4.0,
        'delta': 0.5,
        'theta_uniform': 1.0,
        'Z1': 1,
        'Z2': 1,
        }

        spar = {
        covariance_type = 'ard_se',
        sparse_method = 'uniform',
        n_sparse = 10,
        }

        """

        cur_gap = ''

        if 'name' in desc.keys():
            cur_gap += '%s ' %desc['name']
        else:
            raise ValueError('The name of the descriptor must be set.')

        for key, value in desc.items():
            if key != 'name':
                cur_gap += '%s=%s ' %(key,str(value))

        for key, value in spar.items():
            cur_gap += '%s=%s ' %(key,str(value))

        cur_gap += 'add_species=F'

        self.gaps.append(cur_gap)

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
        command += ' do_copy_at_file=F sparse_separate_file=T gp_file=GAP.xml at_file=train.xyz default_sigma={0.008 0.04 0.04 0} '

        if len(self.gaps) > 0:
            command += 'gap={'
            for i, gap in enumerate(self.gaps):
                command += '%s' %gap
                if i != len(self.gaps)-1:
                    command += ' : '
            command += '}'
        else:
            raise ValueError('There is no gaussian process.')

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
    gf = GapFit(directory='gap-4')

    desc_d = {
            'name': 'distance_2b',
            'cutoff': 4.0,
            'delta': 0.5,
            'theta_uniform': 1.0,
            'Z1': 6,
            'Z2': 8,
    }

    desc_a = {
            'name': 'angle_3b',
            'cutoff': 3,
            'delta': 0.5,
            'theta_uniform': 0.5,
            'Z': 6,
            'Z1': 6,
            'Z2': 8,
    }

    soap = {

            }

    spar = {
            'covariance_type': 'ard_se',
            'sparse_method': 'uniform',
            'n_sparse': 100,
    }

    # =====
    # distance, C-C(x), C-O, C-Pt, O-O, O-Pt, Pt-Pt
    desc_d['cutoff'] = 3.0 # C-O
    gf.add_gap(desc_d,spar)

    desc_d['cutoff'] = 3.0 # C-Pt
    desc_d['Z1'], desc_d['Z2'] = 6, 78
    gf.add_gap(desc_d,spar)

    desc_d['cutoff'] = 2.0 # O-O
    desc_d['Z1'], desc_d['Z2'] = 8, 8
    gf.add_gap(desc_d,spar)

    desc_d['cutoff'] = 3.0 # O-Pt
    desc_d['Z1'], desc_d['Z2'] = 8, 78
    gf.add_gap(desc_d,spar)

    desc_d['cutoff'] = 4.0 # Pt-Pt
    desc_d['Z1'], desc_d['Z2'] = 78, 78
    gf.add_gap(desc_d,spar)

    # auto generation
    elements = [6,8,78]
    nelements = len(elements)

    two_body = []
    for i in range(nelements):
        for j in range(1,nelements):
            two_body.append((i,j))

    three_body = []

    # =====
    # angle: n*n+n*(n*(n-1))/2
    spar['n_sparse'] = 250

    #desc_a['Z'], desc_a['Z1'], desc_a['Z2'] = 6, 6, 6 # C_C-C(x)

    #desc_a['Z'], desc_a['Z1'], desc_a['Z2'] = 6, 6, 8 # C_C-O(x)

    #desc_a['Z'], desc_a['Z1'], desc_a['Z2'] = 6, 6, 78 # C_C-Pt(x)

    desc_a['Z'], desc_a['Z1'], desc_a['Z2'] = 6, 8, 8 # C_O-O
    gf.add_gap(desc_a,spar)

    desc_a['Z'], desc_a['Z1'], desc_a['Z2'] = 6, 8, 78 # C_O-Pt
    gf.add_gap(desc_a,spar)

    desc_a['Z'], desc_a['Z1'], desc_a['Z2'] = 6, 78, 78 # C_Pt-Pt
    gf.add_gap(desc_a,spar)

    #desc_a['Z'], desc_a['Z1'], desc_a['Z2'] = 8, 6, 6 # O-C-C

    desc_a['Z'], desc_a['Z1'], desc_a['Z2'] = 8, 6, 8 # O_C-O
    gf.add_gap(desc_a,spar)

    desc_a['Z'], desc_a['Z1'], desc_a['Z2'] = 8, 6, 78 # O_C-Pt
    gf.add_gap(desc_a,spar)

    #desc_a['Z'], desc_a['Z1'], desc_a['Z2'] = 8, 8, 8 # O_O-O

    #desc_a['Z'], desc_a['Z1'], desc_a['Z2'] = 8, 8, 78 # O_O-Pt

    desc_a['Z'], desc_a['Z1'], desc_a['Z2'] = 8, 78, 78 # O_Pt-Pt
    gf.add_gap(desc_a,spar)

    #desc_a['Z'], desc_a['Z1'], desc_a['Z2'] = 78, 6, 6 # Pt_C-C

    desc_a['Z'], desc_a['Z1'], desc_a['Z2'] = 78, 6, 8 # Pt_C-O
    gf.add_gap(desc_a,spar)

    desc_a['Z'], desc_a['Z1'], desc_a['Z2'] = 78, 6, 78 # Pt_C-Pt
    gf.add_gap(desc_a,spar)

    #desc_a['Z'], desc_a['Z1'], desc_a['Z2'] = 78, 8, 8 # Pt_O-O

    desc_a['Z'], desc_a['Z1'], desc_a['Z2'] = 78, 8, 78 # Pt_O-Pt
    gf.add_gap(desc_a,spar)

    desc_a['Z'], desc_a['Z1'], desc_a['Z2'] = 78, 78, 78 # Pt-Pt-Pt
    gf.add_gap(desc_a,spar)

    # para names
    gf.add_para_names(energy='free_energy',force='forces',virial='virial')

    print(gf.generate_command())
    #gf.run()

