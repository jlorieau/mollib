"""
Test the geometry functions
"""
import unittest

from mollib.core import Molecule, measure_distances, measure_angles


class TestGeometry(unittest.TestCase):
    "Tests the geometry methods."

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

