"""MolLib command line interface
"""
# Author: Justin L Lorieau
# Copyright 2016

from pprint import pprint
from mollib.core import Molecule, measure_angle
from mollib.hydrogenate import add_h
from mollib.hbonds import (find_amide_hbond_partners,
                           find_aliphatic_hbond_partners)


if __name__ == "__main__":
    mol = Molecule('2PTN')
    add_h(mol)
#    mol.write_pdb('2MJB_H.pdb')
    hbonds = find_amide_hbond_partners(mol)
    pprint(hbonds)
    hbonds = find_aliphatic_hbond_partners(mol)
    pprint(hbonds)
    mol.write_pdb('output/2PTN_h.pdb')
