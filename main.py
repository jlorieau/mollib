"""MolLib command line interface
"""
# Author: Justin L Lorieau
# Copyright 2016

from pprint import pprint
import logging
import sys

from mollib.core import Molecule, measure_angle
from mollib.plugins.hydrogenate import add_h
from mollib.plugins.hbonds import (find_amide_hbond_partners,
                                   find_aliphatic_hbond_partners)


def do():
    #logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
    mol = Molecule('2PTN')
    add_h(mol)
    #    mol.write_pdb('2MJB_H.pdb')
    hbonds = find_amide_hbond_partners(mol)
    hbonds = find_aliphatic_hbond_partners(mol)
    pprint(hbonds)
    mol.write_pdb('output/2PTN_h.pdb')

if __name__ == "__main__":
    do()
