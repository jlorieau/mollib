"""
Tests the classification functions
"""

import unittest
import os

from mollib.core import Molecule, load_settings
import mollib.core.settings
from mollib.hbonds import classify_residues, find_hbond_partners, settings


def convert_dssp(string):
    """Convert a DSSP classification (like) string into a classification
    dict that can be used in tests.
    
    Parameters
    ----------
    string: str
        The DSSP (like) string.
    """
    answer_key = {' ': '',
                  'E': settings.major_beta,
                  'H': settings.major_alpha,
                  'G': settings.major_310,
                  'I': settings.major_pi,
                  'a': settings.major_beta_turnI,
                  'b': settings.major_beta_turnIp,
                  'c': settings.major_beta_turnII,
                  'd': settings.major_beta_turnIIp,
                  }

    return {count: answer_key[k] for count, k in enumerate(string, 1)}


class TestHbondClassify(unittest.TestCase):

    def setUp(self):
        """Load the settings."""
        load_settings()

    def test_classify_residues(self):
        "Tests the residue classification from hbonds for specific molecule"


        answer_key = {
            # Hemagglutinin fusion peptide
            '2KXA': convert_dssp(' HHHHHHHHHHH HHHHHHHHH  '),

            # Ubiquitin crystal structure
            '1UBQ': convert_dssp(' EEEEEEaa EEEEEE  aa  HHHHHHHHHHHH   aaa'
                                 'EEEEEbb EEEaa  GGGG   ccEEEEEEE     ')
                      }

        msg = "Residue {} is assigned as '{}', but the test has '{}' assigned"

        for identifier, class_dict in answer_key.items():
            mol = Molecule(identifier)
            classify_residues(mol)

            for residue in mol.residues:
                # Skip HETATM molecules
                if '*' in residue.chain.id:
                    continue

                test_classification = class_dict[residue.number]
                actual_classification = residue.classification[0]

                residue_msg = msg.format(residue, actual_classification,
                                         test_classification)

                self.assertEqual(actual_classification, test_classification,
                                 msg=residue_msg)

    def test_energy_hbond(self):
        """Test the assignment of hydrogen bond energies."""
        # Load the molecule
        mol = Molecule('2KXA')

        # Find hydrogen bonds and classify them
        hbonds = find_hbond_partners(mol)

        # Assert that all of the hbonds have an energy associated and that these
        # are float numbers
        for hbond in hbonds:
            msg = ("The hydrogen bond '{}' does not have an 'energy_hbond' " 
                   "assignment.")
            self.assertTrue(hasattr(hbond, 'energy_hbond'),
                            msg=msg.format(hbond.short_repr()))
            self.assertTrue(isinstance(hbond.energy_hbond, float))

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
            msg = ("Residue '{}' does not have an 'energy_ramachandran' "
                   "assignment")
            self.assertTrue(hasattr(residue, 'energy_ramachandran'),
                            msg=msg.format(residue))
            self.assertIsInstance(residue.energy_ramachandran, float)
