"""
Test functions for fixers.
"""
import unittest
import os

from mollib import Molecule
from mollib.pa.fixers import SignFixer
from mollib.pa.data_readers import read_pa_file


class TestFixers(unittest.TestCase):

    def test_signfixer(self):
        # Load the molecule
        mol = Molecule('2KXA')

        # Load the sign-inverted data
        path = os.path.dirname(os.path.abspath(__file__))
        data_ref = read_pa_file(path + '/data/2kxa_sag.inp')

        # Setup the fixer
        fixer = SignFixer(mol)

        # The data needs to have the N-H RDCs flipped in sign
        data_fixed, fixes = fixer.fix(data_ref)
        self.assertEqual(len(fixes), 1)

        # Now test a dataset that does not need sign flipping
        mol = Molecule('2MJB')

        # Load the sign-inverted data
        path = os.path.dirname(os.path.abspath(__file__))
        data_ref = read_pa_file(path + '/data/ubq_bicelle_hn-c.pa')

        # Setup the fixer
        fixer = SignFixer(mol)

        # The data does not need to be fixed.
        data_fixed, fixes = fixer.fix(data_ref)
        self.assertEqual(len(fixes), 0)
