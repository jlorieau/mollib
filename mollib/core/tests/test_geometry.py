"""
Test the geometry functions
"""
import unittest

from nose.plugins.attrib import attr
import numpy as np

from mollib.core import (Molecule, measure_distances, calc_vector,
                         vector_length, within_distance)


class TestGeometry(unittest.TestCase):
    "Tests the geometry methods."

    def test_calc_vector(self):
        mol = Molecule('2KXA')

        vec = calc_vector(mol['A'][3]['N'].pos, mol['A'][3]['H'].pos,
                          normalize=True)
        length = round(vector_length(vec), 2)
        self.assertEqual([round(i, 2) for i in vec],
                         [-0.39, 0.63, -0.67])
        self.assertEqual(length, 1.00)

        vec = calc_vector(mol['A'][3]['N'].pos, mol['A'][3]['H'].pos,
                          normalize=False)
        length = round(vector_length(vec), 2)
        self.assertEqual([round(i, 2) for i in vec],
                         [-0.38, 0.61, -0.66])
        self.assertEqual(length, 0.98)

        # Test iterables. These aren't supported
        with self.assertRaises(TypeError):
            vec = calc_vector((1.0, 2.0), (0.0, -1.0))

    @attr('bench')
    def test_speed_calc_vector(self):
        """"Tests the speed of the calc_vector function.

        The number of vectors is tuned so that this test takes 1.0s.
        """
        np.random.seed(0)
        for i in range(55000):
            pt1 = np.random.rand(3) * 10. - 5.0
            pt2 = np.random.rand(3) * 10. - 5.0
            vec = calc_vector(pt1, pt2)
            self.assertEqual(len(vec), len(pt1))

    @attr('bench')
    def test_speed_vector_length(self):
        """"Tests the speed of the vector_length function.

        The number of vectors is tuned so that this test takes 1.0s.
        """
        np.random.seed(0)
        for i in range(100000):
            vec1 = np.random.rand(3) * 10. - 5.0
            length = vector_length(vec1)
            self.assertGreaterEqual(length, -8.67)
            self.assertLessEqual(length, 8.67)

    def test_within_distance(self):
        mol = Molecule('2KXA')

        # All atoms within 3A of A5-CA
        atoms = within_distance(mol['A'][5]['CA'], cutoff=3.0)
        self.assertEqual(len(atoms), 13)

        # Exclude hydrogens
        atoms = within_distance(mol['A'][5]['CA'], cutoff=3.0,
                                elements='C|N|O')
        self.assertEqual(len(atoms), 7)

        # Exclude intraresidue
        atoms = within_distance(mol['A'][5]['CA'], cutoff=3.0,
                                elements='C|N|O', exclude_intraresidue=True)
        self.assertEqual(len(atoms), 3)


    @attr('bench')
    def test_speed_within_distance(self):
        """"Tests the speed of the vector_length function.

        The number of vectors is tuned so that this test takes 1.0s.
        """
        mol = Molecule('2KXA')

        for i in range(6):
            for atom in mol.atoms:
                within_distance(atom, cutoff=3.0)

    def test_measure_distances(self):
        mol = Molecule('2MUV')

        # Measure simple distances in a single chain
        dists = measure_distances(mol, '23-N', '25-N', )
        self.assertEqual(len(dists), 1)
        self.assertEqual([(i,j) for i,j,dist in dists],
                         [(mol['A'][23]['N'], mol['A'][25]['N'])])

        # Measure multiple distances in a single chain
        dists = measure_distances(mol, '23:25-N', '23:25-N')
        self.assertEqual(len(dists), 3)
        self.assertEqual([(i, j) for i, j, dist in dists],
                         [(mol['A'][23]['N'], mol['A'][24]['N']),
                          (mol['A'][23]['N'], mol['A'][25]['N']),
                          (mol['A'][24]['N'], mol['A'][25]['N']), ])

        # Measure specific residue number deltas
        dists = measure_distances(mol, '23:25-N', '23:25-N', residue_delta=2)
        self.assertEqual(len(dists), 1)
        self.assertEqual([(i, j) for i, j, dist in dists],
                         [(mol['A'][23]['N'], mol['A'][25]['N']), ])

        # Try different atoms
        dists = measure_distances(mol, '23:25-N', '23:25-CA')
        self.assertEqual(len(dists), 9)
        self.assertEqual([(i, j) for i, j, dist in dists],
                         [(mol['A'][23]['N'], mol['A'][23]['CA']),
                          (mol['A'][23]['N'], mol['A'][24]['CA']),
                          (mol['A'][23]['N'], mol['A'][25]['CA']),
                          (mol['A'][24]['N'], mol['A'][23]['CA']),
                          (mol['A'][24]['N'], mol['A'][24]['CA']),
                          (mol['A'][24]['N'], mol['A'][25]['CA']),
                          (mol['A'][25]['N'], mol['A'][23]['CA']),
                          (mol['A'][25]['N'], mol['A'][24]['CA']),
                          (mol['A'][25]['N'], mol['A'][25]['CA']), ])

        # Exclude intra residue items
        dists = measure_distances(mol, '23:25-N', '23:25-CA',
                                  exclude_intra=True)
        self.assertEqual(len(dists), 6)
        self.assertEqual([(i, j) for i, j, dist in dists],
                         [(mol['A'][23]['N'], mol['A'][24]['CA']),
                          (mol['A'][23]['N'], mol['A'][25]['CA']),
                          (mol['A'][24]['N'], mol['A'][23]['CA']),
                          (mol['A'][24]['N'], mol['A'][25]['CA']),
                          (mol['A'][25]['N'], mol['A'][23]['CA']),
                          (mol['A'][25]['N'], mol['A'][24]['CA']), ])

        # Try measurements from different chains. The same atoms (within a
        # chain) are removed
        dists = measure_distances(mol, 'A:D.23-CA', 'A:D.23-CA')
        self.assertEqual(len(dists), 6)
        self.assertEqual([(i, j) for i, j, dist in dists],
                         [(mol['A'][23]['CA'], mol['B'][23]['CA']),
                          (mol['A'][23]['CA'], mol['C'][23]['CA']),
                          (mol['A'][23]['CA'], mol['D'][23]['CA']),
                          (mol['B'][23]['CA'], mol['C'][23]['CA']),
                          (mol['B'][23]['CA'], mol['D'][23]['CA']),
                          (mol['C'][23]['CA'], mol['D'][23]['CA']), ])

        # Only include measurements within the same chain. This returns nothing
        # since all measurements are between chains
        dists = measure_distances(mol, 'A:D.23-CA', 'A:D.23-CA',
                                  only_intra_chain=True)
        self.assertEqual(len(dists), 0)

        # Test chain and residue ranges, only show residue delta of 1
        dists = measure_distances(mol, 'A:D.23:24-CA', 'A:D.23:24-CA',
                                  residue_delta=1)
        self.assertEqual(len(dists), 16)
        self.assertEqual([(i, j) for i, j, dist in dists],
                         [(mol['A'][23]['CA'], mol['A'][24]['CA']),
                          (mol['A'][23]['CA'], mol['B'][24]['CA']),
                          (mol['A'][23]['CA'], mol['C'][24]['CA']),
                          (mol['A'][23]['CA'], mol['D'][24]['CA']),
                          (mol['B'][23]['CA'], mol['A'][24]['CA']),
                          (mol['B'][23]['CA'], mol['B'][24]['CA']),
                          (mol['B'][23]['CA'], mol['C'][24]['CA']),
                          (mol['B'][23]['CA'], mol['D'][24]['CA']),
                          (mol['C'][23]['CA'], mol['A'][24]['CA']),
                          (mol['C'][23]['CA'], mol['B'][24]['CA']),
                          (mol['C'][23]['CA'], mol['C'][24]['CA']),
                          (mol['C'][23]['CA'], mol['D'][24]['CA']),
                          (mol['D'][23]['CA'], mol['A'][24]['CA']),
                          (mol['D'][23]['CA'], mol['B'][24]['CA']),
                          (mol['D'][23]['CA'], mol['C'][24]['CA']),
                          (mol['D'][23]['CA'], mol['D'][24]['CA']),
                          ])

