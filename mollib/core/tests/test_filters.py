import unittest
from mollib.core import Molecule, filter_atoms


class TestAtomFilters(unittest.TestCase):
    "Tests the atom filter functions."

    def test_filter_atoms(self):
        "Test the utils.filter_atoms function."
        # A protein structure with multiple chains
        mol = Molecule('2MUV')

        # Test the exclude_intra filter
        self.assertFalse(filter_atoms(mol['A'][23]['N'],
                                      mol['A'][23]['CA'],
                                      mol['A'][23]['C']))

        self.assertTrue(filter_atoms(mol['A'][23]['N'],
                                     mol['A'][23]['CA'],
                                     mol['A'][23]['C'], exclude_intra=True))

        self.assertFalse(filter_atoms(mol['A'][23]['N'],
                                      mol['A'][23]['CA'],
                                      mol['A'][24]['C'], exclude_intra=True))

        self.assertFalse(filter_atoms(mol['A'][23]['N'],
                                      mol['A'][23]['CA'],
                                      mol['B'][23]['C'], exclude_intra=True))

        # Test the only_intra filter
        self.assertFalse(filter_atoms(mol['A'][23]['N'],
                                      mol['A'][23]['CA'],
                                      mol['A'][23]['C'], only_intra=True))

        self.assertTrue(filter_atoms(mol['A'][23]['N'],
                                      mol['A'][23]['CA'],
                                      mol['A'][24]['C'], only_intra=True))

        self.assertTrue(filter_atoms(mol['A'][23]['N'],
                                     mol['A'][23]['CA'],
                                     mol['B'][23]['C'], only_intra=True))

        # Test the only_intra_chain filter
        self.assertFalse(filter_atoms(mol['A'][23]['N'],
                                      mol['A'][23]['CA'],
                                      mol['A'][23]['C'], only_intra_chain=True))

        self.assertFalse(filter_atoms(mol['A'][23]['N'],
                                     mol['A'][23]['CA'],
                                     mol['A'][24]['C'], only_intra_chain=True))

        self.assertTrue(filter_atoms(mol['A'][23]['N'],
                                     mol['A'][23]['CA'],
                                     mol['B'][23]['C'], only_intra_chain=True))

        # Test the exclude_intra_chain filter
        self.assertTrue(filter_atoms(mol['A'][23]['N'],
                                      mol['A'][23]['CA'],
                                      mol['A'][23]['C'],
                                      exclude_intra_chain=True))

        self.assertTrue(filter_atoms(mol['A'][23]['N'],
                                      mol['A'][23]['CA'],
                                      mol['A'][24]['C'],
                                      exclude_intra_chain=True))

        self.assertFalse(filter_atoms(mol['A'][23]['N'],
                                     mol['A'][23]['CA'],
                                     mol['B'][23]['C'],
                                     exclude_intra_chain=True))

        # Test the residue_delta filter
        self.assertTrue(filter_atoms(mol['A'][23]['N'],
                                     mol['A'][23]['CA'],
                                     mol['A'][23]['C'],
                                     residue_delta=1))

        self.assertFalse(filter_atoms(mol['A'][23]['N'],
                                      mol['A'][23]['CA'],
                                      mol['A'][24]['C'],
                                      residue_delta=1))

        self.assertTrue(filter_atoms(mol['A'][23]['N'],
                                     mol['A'][23]['CA'],
                                     mol['B'][23]['C'],
                                     residue_delta=1))

        # This filter ignores chains
        self.assertFalse(filter_atoms(mol['A'][23]['N'],
                                      mol['A'][23]['CA'],
                                      mol['B'][24]['C'],
                                      residue_delta=1))

        # Test the bonded filter
        self.assertFalse(filter_atoms(mol['A'][23]['N'],
                                     mol['A'][23]['CA'],
                                     mol['A'][23]['C'],
                                     bonded=True))

        self.assertTrue(filter_atoms(mol['A'][23]['N'],
                                     mol['A'][23]['CA'],
                                     mol['A'][24]['C'],
                                     bonded=True))

        self.assertTrue(filter_atoms(mol['A'][23]['N'],
                                      mol['A'][23]['CA'],
                                      mol['B'][23]['C'],
                                      bonded=True))