import unittest
from mollib.core import Molecule, measure_angle
from mollib.hydrogenate import add_h


def in_range(value, target, tolerance):
    return (target - tolerance <= value <= target + tolerance)


class TestMolLib(unittest.TestCase):
    """Tests hydrogenation functions"""

    def setUp(self):
        mol = Molecule('2KXA')
        mol.strip_atoms(element='H')
        add_h(mol)
        self.mol = mol

    def test_add_two_sp3_h(self):
        """Tests the add_two_sp3_h function."""
        tolerance = 0.3  # degrees

        # Measure the HA2-CA-HA3 angle to make sure they're all tetrahedral
        for res in self.mol.residues:
            if res.name != 'GLY':
                continue
            angle = measure_angle(res['HA2'], res['CA'], res['HA3'])
            self.assertTrue(in_range(angle, 109.5, tolerance))
