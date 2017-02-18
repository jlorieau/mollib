import unittest

from mollib import Molecule
from mollib.pa.process_molecule import ProcessNHDipole


class TestProcess(unittest.TestCase):
    """Test the processing of molecules in dipoles and anistropic chemical
    shifts"""

    def test_NH_dipole(self):
        """Tests the construction of NH dipoles."""

        mol = Molecule('2MJB')
        process = ProcessNHDipole(mol)
        dipole_arrays = process.process()

        # The dipole arrays should be a list of dicts. There should be one
        # item in the list because we specified only one molecule.
        self.assertEqual(1, len(dipole_arrays))

        # The number of items in the dipolar array dict should be equal to the
        # number of backbone NH dipoles in the molecule.
        no_NHs = sum([1 if 'H' in res else 0 for res in mol.residues])
        self.assertEqual(len(dipole_arrays[0]), no_NHs)

        # All of the arrays should contain five elements
        for arr in dipole_arrays[0].values():
            self.assertEqual(arr.size, 5)


