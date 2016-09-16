"""
Tests for hbonds
"""
import unittest

from nose.plugins.attrib import attr

from mollib import Molecule
from mollib.hbonds import find_hbond_partners


class TestHydrogenBonds(unittest.TestCase):
    "Tests the objects and functions for hbonds."

    @attr('bench')
    def test_large_molecule(self):
        """Tests the detection of hbonds in a large protein.

        .. note:: This test takes about 38.43s. The cProfile bottlenecks are:
                  - 16.134s : mollib.core.vector_length
                  - 8.200s  : mollib.hbonds.dipole_distances
                  - 7.663s  : mollib.core.calc_vector
                  - 2.861s  : mollib.core.measure_distance

                  The following optimizations had the following test times:
                  - 23.4s   : convert vector_length, calc_vector to cython
                  - 21.2s   : optimize vector_length in cython
                  -  5.1s   : cythonize measure_distance and cut down the
                              number of call to it.
        """

        mol = Molecule('2N18')
        hbonds = find_hbond_partners(mol)

        self.assertEqual(len(hbonds), 378)

    
