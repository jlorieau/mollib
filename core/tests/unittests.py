"""Unittests for the core module.
"""
# Author: Justin L Lorieau
# Copyright 2016

import unittest
from datetime import datetime

from mollib.core import Molecule

class TestMolLib(unittest.TestCase):

    performance_tests = False

    def test_large_molecule(self):
        "Tests the parsing and performance of a very large protein complex."
        import string

        if self.performance_tests:
            # 5.5s - default atom generator
            # 6.0s - swithching to atom generator 2
            # 3.6s - after switching convert function to regexes. (2.75MB/s)
            id = '3H0G'  # RNA Polymerase II from Schizosaccharomyces pombe
            start = datetime.now()
            Molecule(id)
            stop = datetime.now()
            time = (stop - start).total_seconds()
            print("Loaded {id} in {time} seconds".format(id=id, time=time))

        mol = Molecule('3H0G')

        # Test that all of the chains were read in correctly
        chains = list(string.ascii_uppercase)[:24]  # chains A-X
        chains += ['A*', 'B*', 'C*', 'I*', 'J*', 'L*', 'M*', 'N*', 'O*', 'U*',
                   'V*', 'X*']
        self.assertEqual([c.id for c in mol.chains], sorted(chains))
        self.assertEqual(mol.chain_size, len(chains))

        # Test the molecular mass of each chain
        for chain in mol.chains:
            self.assertGreater(chain.mass, 0.)

        self.assertAlmostEqual(mol.mass, 833388.28, 2)

    def test_multiple_models(self):
        """Tests reading PDB files with multiple models. Only the first model
        should be read in."""
        mol = Molecule('2KXA')  # 20 models

        # These are the coordinates for this atom of the first model
        self.assertEqual(mol['A'][3]['N'].pos[0], 13.766)
        self.assertEqual(mol['A'][3]['N'].pos[1], -3.965)
        self.assertEqual(mol['A'][3]['N'].pos[2], 5.893)

    def test_residue_ordering(self):
        """Tests the linked lists of residues."""
        mol = Molecule('2KXA')

        last_residues = [r.last_residue.number
                         if r.last_residue is not None else None
                         for r in mol.residues]
        self.assertEqual(last_residues,
                         [None, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13,
                          14, 15, 16, 17, 18, 19, 20, 21, 22, 23])
