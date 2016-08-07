import unittest
from mollib.core import Molecule, measure_angle, measure_dihedral
from mollib.hydrogenate import add_h, methines, methylenes


def in_range(value, target, tolerance, print_false=True):
    r = (target - tolerance <= value <= target + tolerance)
    if r is False and print_false:
        print('Value {} is not within {} +/- {}'.format(value, target,
                                                        tolerance))
    return r


class TestMolLib(unittest.TestCase):
    """Tests hydrogenation functions"""

    def setUp(self):
        mol = Molecule('2PTN')
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

    def test_methines(self):
        """Tests the geometry of methines added."""
        tolerance = 11. # degrees

        residue_found = {k:False for k in methylenes}
        for residue in self.mol.residues:
            if residue.name in methines:
                for atom_name, target_name, a1_name, a2_name, a3_name in \
                    methines[residue.name]:

                    h = residue[atom_name]
                    target_atom = residue[target_name]
                    a1 = residue[a1_name]
                    a2 = residue[a2_name]
                    a3 = residue[a3_name]

                    for a in (a1, a2, a3):
                        angle = measure_angle(h, target_atom, a)
                        msg = "{}_{}_{} ({}) not within {} +/- {}"
                        msg = msg.format(h, target_atom, a, angle,
                                         109.5, tolerance)
                        self.assertTrue(in_range(angle, 109.5, tolerance), msg)

                    # This type of residue has been found
                    residue_found[residue.name] = True
        self.assertTrue(all(residue_found))

    def test_methylenes(self):
        """Tests the geometry and stereospecificity of methylene hydrogens
        added."""
        tolerance = 10.  # degrees

        residue_found = {k: False for k in methylenes}
        for residue in self.mol.residues:
            if residue.name in methylenes:
                for atom_name, target_name, a1_name, a2_name in \
                        methylenes[residue.name]:

                    h2 = residue[atom_name + '2']
                    h3 = residue[atom_name + '3']

                    target_atom = residue[target_name]
                    a1 = residue[a1_name]
                    a2 = residue[a2_name]

                    # Test that the new methylene is tetrahedral
                    for h,a in ((h2, a1), (h2, a2), (h3, a1), (h3, a2)):
                        angle = measure_angle(h, target_atom, a)
                        msg = "{}_{}_{} ({}) not within {} +/- {}"
                        msg = msg.format(h, target_atom, a, angle,
                                         109.5, tolerance)
                        self.assertTrue(in_range(angle, 109.5, tolerance),
                                        msg)

                    # Test that the assignment respects the stereospecificity
                    # H2 is Pro-R and H3 is Pro-S
                    proR_angle = measure_dihedral(target_atom, a1, a2, h2)
                    proS_angle = measure_dihedral(target_atom, a1, a2, h3)

                    msg = "{}_{} ({}) stereospecific assignment error."

                    self.assertTrue(in_range(proR_angle, -30., tolerance),
                                    msg.format(h2, target_atom, proR_angle))
                    self.assertTrue(in_range(proS_angle, 30., tolerance),
                                    msg.format(h3, target_atom, proS_angle))

                    # This type of residue has been found
                    residue_found[residue.name] = True
        self.assertTrue(all(residue_found.values()))
