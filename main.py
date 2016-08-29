"""MolLib command line interface
"""
# Author: Justin L Lorieau
# Copyright 2016

from pprint import pprint
import logging
import sys

from mollib.core import Molecule, measure_angle
from mollib.hydrogens import add_hydrogens, add_one_sp2_h
#from mollib.hydrogenate import add_h
from mollib.hbonds import (find_amide_hbond_partners,
                           find_aliphatic_hbond_partners)


def do():
    #logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
    #mol = Molecule('2MUV')
    mol = Molecule('2PTN')
    mol = Molecule('2MJB')
    add_hydrogens(mol)
    #    mol.write_pdb('2MJB_H.pdb')
    #hbonds = find_amide_hbond_partners(mol)
    #hbonds = find_aliphatic_hbond_partners(mol)
    #pprint(hbonds)
    #mol.write_pdb('output/2PTN_h.pdb')
    #mol.write_pdb('output/2MUV_h.pdb')
    mol.write_pdb('output/2MJB_h.pdb')

if __name__ == "__main__":
    do()
