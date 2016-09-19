import unittest

from mollib import Molecule
from mollib.hydrogens import add_hydrogens


class TestHydrogenateMolecules(unittest.TestCase):
    "Test the hydrogenation of specific molecules"

    def test_1C1D(self):
        mol = Molecule('1C1D')
        add_hydrogens(mol)