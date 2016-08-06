"""MolLib command line interface

   @Author:             Justin L Lorieau <jlorieau>
   @Date:               2016-08-06T06:07:52-05:00
   @Last modified by:   jlorieau
   @Last modified time: 2016-08-06T06:09:30-05:00
   @License:            Copyright 2016
"""
from pprint import pprint
from mollib import Molecule
from mollib.protonate import add_h
from mollib.hbonds import (find_amide_hbond_partners,
                           find_aliphatic_hbond_partners)


if __name__ == "__main__":
    mol = Molecule('2N7J')
    add_h(mol)
#    mol.write_pdb('2MJB_H.pdb')
    hbonds = find_amide_hbond_partners(mol)
    pprint(hbonds)
    hbonds = find_aliphatic_hbond_partners(mol)
    pprint(hbonds)
