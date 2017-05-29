"""
Test the geometry functions
"""
import unittest
from math import sqrt

import numpy as np

from mollib.core import (Molecule, measure_distances, calc_vector,
                         vector_length, within_distance)
from mollib.core.geometry_box import Box


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

    def test_measure_distances(self):
        mol = Molecule('2MUV')

        # Measure simple distances in a single chain
        dists = measure_distances(mol, '23.N', '25.N', )
        self.assertEqual(len(dists), 1)
        self.assertEqual([(i,j) for i,j,dist in dists],
                         [(mol['A'][23]['N'], mol['A'][25]['N'])])

        # Measure multiple distances in a single chain
        dists = measure_distances(mol, '23:25.N', '23:25.N')
        self.assertEqual(len(dists), 3)
        self.assertEqual([(i, j) for i, j, dist in dists],
                         [(mol['A'][23]['N'], mol['A'][24]['N']),
                          (mol['A'][23]['N'], mol['A'][25]['N']),
                          (mol['A'][24]['N'], mol['A'][25]['N']), ])

        # Measure specific residue number deltas
        dists = measure_distances(mol, '23:25.N', '23:25.N', residue_delta=2)
        self.assertEqual(len(dists), 1)
        self.assertEqual([(i, j) for i, j, dist in dists],
                         [(mol['A'][23]['N'], mol['A'][25]['N']), ])

        # Try different atoms
        dists = measure_distances(mol, '23:25.N', '23:25.CA')
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
        dists = measure_distances(mol, '23:25.N', '23:25.CA',
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
        dists = measure_distances(mol, 'A:D.23.CA', 'A:D.23.CA')
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
        dists = measure_distances(mol, 'A:D.23.CA', 'A:D.23.CA',
                                  only_intra_chain=True)
        self.assertEqual(len(dists), 0)

        # Test chain and residue ranges, only show residue delta of 1
        dists = measure_distances(mol, 'A:D.23:24.CA', 'A:D.23:24.CA',
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

    def test_box(self):
        """Tests the Box class for storing and accessing nearest-neighbor
        points."""
        # create the points
        point_range = (-20., 20.)
        no_points = 5000
        node_size = 2.0

        np.random.seed(0)
        points = []
        for i in range(no_points):
            points.append(np.random.rand(3)*point_range[1] - point_range[1]/2.)

        # Setup the distance function
        def dist(p1, p2):
            return sqrt(sum([(i-j)**2 for i,j in zip(p1, p2)]))

        # Create the box
        box = Box(points, node_size=node_size)
        self.assertEqual(len(points), box.num_points())

        # Test that the points are correctly retrieved and compare to a brute
        # search
        closest_points = list(box.get_points(points[0], radius=2.0))
        brute_closest = [p for p in points if dist(p, points[0]) < 2.0]

        self.assertGreater(len(closest_points), 0)
        self.assertEqual(len(closest_points), len(brute_closest))

        # Access an existing node
        node = box.get_node(points[0])
        self.assertIn(points[0], node)

        # Access an non-existing node. This raises an IndexError
        fake_point = (point_range[0]*3., point_range[0]*3., point_range[0]*3)
        with self.assertRaises(IndexError):
            node = box.get_node(fake_point)

        # Try making a box with mismatched point dimensions. This raises an
        # index error
        points2 = [(1., 2., 3.),
                   (4., 5.,)]
        with self.assertRaises(IndexError):
            box = Box(points2)

        # Try making a box with a negative node_size. This raises an assertion
        # error
        with self.assertRaises(AssertionError):
            box = Box(points, node_size=-2.)