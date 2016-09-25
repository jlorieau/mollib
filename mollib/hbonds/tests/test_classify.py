"""
Tests the classification functions
"""

import unittest

from mollib import Molecule
from mollib.hbonds import classify_residues


class TestHbondClassify(unittest.TestCase):

    def test_classify_residues(self):
        "Tests the residue classification from hbonds for specific molecule"

        # Hemagglutinin fusion peptide
        answer_key = {'2KXA': {1:  'isolated',
                               2:  'alpha-helix',
                               3:  'alpha-helix',
                               4:  'alpha-helix',
                               5:  'alpha-helix',
                               6:  'alpha-helix',
                               7:  'alpha-helix',
                               8:  'alpha-helix',
                               9:  'alpha-helix',
                               10: 'alpha-helix',
                               11: 'alpha-helix',
                               12: '',
                               13: 'isolated',
                               14: 'alpha-helix',
                               15: 'alpha-helix',
                               16: 'alpha-helix',
                               17: 'alpha-helix',
                               18: 'alpha-helix',
                               19: 'alpha-helix',
                               20: 'alpha-helix',
                               21: 'alpha-helix',
                               22: 'alpha-helix',
                               23: '',
                               24: 'isolated',
                               }
                      }

        for identifier, class_dict in answer_key.items():
            mol = Molecule(identifier)
            classify_residues(mol)

            for residue in mol.residues:
                res_class = class_dict[residue.number]
                self.assertEqual(residue.hbond_classification,
                                 res_class)

