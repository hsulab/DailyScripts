#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

def calc_water_num(volume):
    """ Calculate number of water molecules at 298.15K with given volume
        (AA^3)
    """
    # water density
    water_molecule_weight = 18.0152 # g/mol
    water_density = 0.997074 # g/cm^3 298.15K
    n_avogadro = 6.02e23 # mol
    # volume = 16*13*10 # Å^3

    n_water_molecule_per_A = (water_density / water_molecule_weight) * n_avogadro * 1e-24

    return np.floor(n_water_molecule_per_A * volume)

print(calc_water_num(12.5**3))

if __name__ == "__main__":
    pass
