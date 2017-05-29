import unittest

from mollib import Molecule
from mollib.pa.process_molecule import Process


class TestProcess(unittest.TestCase):
    """Test the processing of molecules in dipoles and anistropic chemical
    shifts"""

    def test_dipole(self):
        """Tests the construction of NH dipoles."""

        mol = Molecule('2MJB')
        process = Process(mol)
        interactions = process.process()

        # The interactions should be a list of dicts. There should be one
        # item in the list because we specified only one molecule.
        self.assertEqual(1, len(interactions))

        # The number of items in the interactions dict should be equal to the
        # of interactions specified (more than one)
        self.assertGreater(len(interactions[0]), 1)

        # All of the interaction arrays should contain five elements
        for scale, arr in interactions[0].values():
            self.assertEqual(arr.size, 5)


