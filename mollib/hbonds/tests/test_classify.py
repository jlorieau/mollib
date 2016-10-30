"""
Tests the classification functions
"""

import unittest
import os

from mollib.core import Molecule, load_settings
import mollib.core.settings
from mollib.hbonds import classify_residues


class TestHbondClassify(unittest.TestCase):

    def setUp(self):
        "Load the settings."
        load_settings()

    def test_classify_residues(self):
        "Tests the residue classification from hbonds for specific molecule"

        # Hemagglutinin fusion peptide
        answer_key = {'2KXA': {1:  ('isolated', ''),
                               2:  ('alpha-helix', 'N-term'),
                               3:  ('alpha-helix', 'N-term'),
                               4:  ('alpha-helix', ''),
                               5:  ('alpha-helix', ''),
                               6:  ('alpha-helix', ''),
                               7:  ('alpha-helix', ''),
                               8:  ('alpha-helix', ''),
                               9:  ('alpha-helix', ''),
                               10: ('alpha-helix', ''),
                               11: ('alpha-helix', 'C-term'),
                               12: ('alpha-helix', 'C-term'),
                               13: ('isolated', ''),
                               14: ('alpha-helix', 'N-term'),
                               15: ('alpha-helix', 'N-term'),
                               16: ('alpha-helix', ''),
                               17: ('alpha-helix', ''),
                               18: ('alpha-helix', ''),
                               19: ('alpha-helix', ''),
                               20: ('alpha-helix', ''),
                               21: ('alpha-helix', 'C-term'),
                               22: ('alpha-helix', 'C-term'),
                               23: ('isolated', ''),
                               24: ('isolated', ''),
                               }
                      }

        for identifier, class_dict in answer_key.items():
            mol = Molecule(identifier)
            classify_residues(mol)

            for residue in mol.residues:
                items = class_dict[residue.number]
                hbond_classification, hbond_modifier = items
                print(residue, residue.hbond_classification, residue.hbond_modifier)
                self.assertEqual(residue.hbond_classification,
                                 hbond_classification)
                self.assertEqual(residue.hbond_modifier,
                                 hbond_modifier)

    def test_energy_ramachandran(self):
        """Test the 'energy_ramachandran` property of residues, set by
        classify_residues.
        """
        # Load the molecule
        mol = Molecule('2KXA')

        # First confuse the path for the ramachandran_dataset_path so that the
        # datasets cannot be found.
        correct_path = mollib.core.settings.ramachandran_dataset_path
        mollib.core.settings.ramachandran_dataset_path = ''

        # Try classify_residues. All of the 'energy_ramachandran' attributes
        # should not be assigned because the datasets could not be found.
        for residue in mol.residues:
            self.assertFalse(hasattr(residue, 'energy_ramachandran'))

        # With the correct path, the datasets should be found and the energies
        # are correctly set
        correct_path = os.path.join('../..', correct_path)
        mollib.core.settings.ramachandran_dataset_path = correct_path

        classify_residues(mol)

        # The 'energy_ramachandran' attributes should now be assigned and
        # have float values
        for residue in mol.residues:
            self.assertTrue(hasattr(residue, 'energy_ramachandran'))
            self.assertIsInstance(residue.energy_ramachandran, float)
