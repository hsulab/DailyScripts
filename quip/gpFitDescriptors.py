#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from ase.data import atomic_numbers

class CovarianceWrapper(object):
    def __init__(self, covariance_type, delta, theta):
        # check
        implemented_covariances = ['ard_se']
        if covariance_type not in implemented_covariances:
            raise ValueError('Covariance %s not implented!' %covariance_type)
        else:
            self.covariance_type = covariance_type

        self.delta = delta
        self.theta = theta

        pass

    def __str__(self):
        command = ('covariance_type={:s} delta={:f} theta_uniform={:f}')\
                .format(self.covariance_type, self.delta, self.theta)
        return command

class SparseWrapper(object):
    def __init__(self, sparse_method, n_sparse):
        implemented_sparses = ['uniform']
        if sparse_method not in implemented_sparses:
            raise ValueError('Sparse %s not implented!' %sparse_method)
        else:
            self.sparse_method = sparse_method
        pass

        self.n_sparse = n_sparse

    def __str__(self):
        return 'sparse_method=%s n_sparse=%d' %(self.sparse_method, self.n_sparse)

class DescriptorWrapper(object):
    def __init__(self,name,cutoff,elements=[]):
        # check name
        implemented_descriptors = ['distance_2b', 'angle_3b']
        if name not in implemented_descriptors:
            raise ValueError('Descriptor %s not implented!' %name)
        else:
            self.name = name

        # check cutoff
        if cutoff < 0:
            raise ValueError('Cutoff must be positive')
        else:
            self._cutoff = cutoff

        # symbols and numbers
        if elements:
            for elem in elements:
                if elem not in atomic_numbers.keys():
                    raise ValueError('No such chemicl symbol %s' %elem)
            else:
                self._chemical_symbols = elements
        else:
            if self.name == 'distance_2b':
                self._chemical_symbols = ['X']*2
            elif self.name == 'angle_3b':
                self._chemical_symbols = ['X']*3
            else:
                pass # already check at first
        self._chemical_numbers = self.symbol2number(self._chemical_symbols)

        pass

    # ===== properties
    @property
    def cutoff(self):
        """Get the current cutoff."""
        return self._cutoff

    @cutoff.setter
    def cutoff(self, value):
        self._cutoff = value
        return

    @property
    def command(self):
        self._command = '%s cutoff=%f' %(self.name, self._cutoff)
        return self._command

    @property
    def chemical_symbols(self):
        """Chemical symbols"""
        return self._chemical_symbols

    @chemical_symbols.setter
    def chemical_symbols(self,symbols):
        # update numbers as well
        self._chemical_symbols = symbols
        self._chemical_numbers = self.symbol2number(symbols)
        return 

    # ===== methods
    @staticmethod
    def symbol2number(symbols):
        numbers = []
        for sym in symbols:
            numbers.append(atomic_numbers[sym])

        return numbers

    def __str__(self):
        command = '%s cutoff=%f ' %(self.name, self._cutoff)
        if self.name == 'distance_2b':
            command += 'Z1={:d} Z2={:d}'.format(*self._chemical_numbers)
        elif self.name == 'angle_3b':
            command += 'Z={:d} Z1={:d} Z2={:d}'.format(*self._chemical_numbers)
        else:
            pass

        return command

def concatenate_descriptor(desc, cov, spar, add_species=False):
    if add_species:
        command = ' '.join([desc.command,str(cov),str(spar)])
        command += ' add_species=T'
    else:
        if 'X' in desc.chemical_symbols:
            raise ValueError('Chemical symbols not set')
        command = ' '.join([str(desc),str(cov),str(spar)])
        command += ' add_species=F'

    return command

# templates
dis_temp = DescriptorWrapper('distance_2b', 4.0, [])
ang_temp = DescriptorWrapper('angle_3b', 4.0, [])

cov = CovarianceWrapper('ard_se', 0.5, 0.5)

spar_050 = SparseWrapper('uniform', 50)
spar_100 = SparseWrapper('uniform', 100)
spar_150 = SparseWrapper('uniform', 150)
spar_250 = SparseWrapper('uniform', 250)
spar_500 = SparseWrapper('uniform', 500)
spar_1000 = SparseWrapper('uniform', 250)

# descriptor list
desc_list = []

distances = [
        ['C','O'],['C','Pt'],['O','O'],['O','Pt'],['Pt','Pt'],
        ]

# distance, C-C(x), C-O, C-Pt, O-O, O-Pt, Pt-Pt
#dis_temp.chemical_symbols = ['C','O']; dis_temp.cutoff = 3.0
#desc_list.append(concatenate_descriptor(dis_temp,cov,spar_500))
#
#dis_temp.chemical_symbols = ['C','Pt']; dis_temp.cutoff = 3.0
#desc_list.append(concatenate_descriptor(dis_temp,cov,spar_500))
#
#dis_temp.chemical_symbols = ['O','O']; dis_temp.cutoff = 3.0
#desc_list.append(concatenate_descriptor(dis_temp,cov,spar_500))
#
#dis_temp.chemical_symbols = ['O','Pt']; dis_temp.cutoff = 3.0
#desc_list.append(concatenate_descriptor(dis_temp,cov,spar_500))
#
#dis_temp.chemical_symbols = ['Pt','Pt']; dis_temp.cutoff = 4.0
#desc_list.append(concatenate_descriptor(dis_temp,cov,spar_500))
for symbols in distances:
    dis_temp.chemical_symbols = symbols
    dis_temp.cutoff = 3.0
    if symbols == ['Pt','Pt']:
        dis_temp.cutoff = 4.0
    desc_list.append(concatenate_descriptor(dis_temp,cov,spar_050))

# angle
# C_C-C(x), C_C-O(x), C_C-Pt(x), C_O-O, C_O-Pt, C_Pt-Pt
# O_C-C(x), O_C-O, O_C-Pt, O_O-O(x), O_O-Pt(x), O_Pt-Pt
# Pt_C-C(x), Pt_C-O, Pt_C-Pt, Pt_O-O(x), Pt_O-Pt, Pt-Pt-Pt
ang_temp.chemical_symbols = ['C','O','O']; ang_temp.cutoff = 3.0
desc_list.append(concatenate_descriptor(ang_temp,cov,spar_1000))

ang_temp.chemical_symbols = ['C','O','Pt']; ang_temp.cutoff = 3.0
desc_list.append(concatenate_descriptor(ang_temp,cov,spar_1000))

ang_temp.chemical_symbols = ['C','Pt','Pt']; ang_temp.cutoff = 3.0
desc_list.append(concatenate_descriptor(ang_temp,cov,spar_1000))

ang_temp.chemical_symbols = ['O','C','O']; ang_temp.cutoff = 3.0
desc_list.append(concatenate_descriptor(ang_temp,cov,spar_1000))

ang_temp.chemical_symbols = ['O','C','Pt']; ang_temp.cutoff = 3.0
desc_list.append(concatenate_descriptor(ang_temp,cov,spar_1000))

ang_temp.chemical_symbols = ['O','Pt','Pt']; ang_temp.cutoff = 3.0
desc_list.append(concatenate_descriptor(ang_temp,cov,spar_1000))

ang_temp.chemical_symbols = ['Pt','C','O']; ang_temp.cutoff = 3.0
desc_list.append(concatenate_descriptor(ang_temp,cov,spar_1000))

ang_temp.chemical_symbols = ['Pt','C','Pt']; ang_temp.cutoff = 3.0
desc_list.append(concatenate_descriptor(ang_temp,cov,spar_1000))

ang_temp.chemical_symbols = ['Pt','O','Pt']; ang_temp.cutoff = 3.0
desc_list.append(concatenate_descriptor(ang_temp,cov,spar_1000))

ang_temp.chemical_symbols = ['Pt','Pt','Pt']; ang_temp.cutoff = 4.0
desc_list.append(concatenate_descriptor(ang_temp,cov,spar_1000))

#print(desc_list)
