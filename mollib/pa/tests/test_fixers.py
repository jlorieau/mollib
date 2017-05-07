"""
Test functions for fixers.
"""
import unittest
import os

from mollib import Molecule
from mollib.pa.fixers import SignFixer, OutlierFixer, NHScaleFixer
from mollib.pa.data_readers import read_pa_file


class TestFixers(unittest.TestCase):

    def test_signfixer(self):
        """Test the SignFixer."""

        # Load the molecule
        mol = Molecule('2KXA')

        # Load the sign-inverted data
        path = os.path.dirname(os.path.abspath(__file__))
        data_ref = read_pa_file(path + '/data/2kxa_dGpG_nh-caha.pa')

        # Setup the fixer
        fixer = SignFixer(mol)

        # The data needs to have the N-H RDCs flipped in sign
        data_fixed, fixes = fixer.fix(data_ref)
        self.assertEqual(len(fixes), 1)

        # Check to see if the overall Q-factor is good
        RMS, Q, data_pred = fixer.fit(data_fixed)

        self.assertLessEqual(Q, 15.0)

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

    def test_nhscalefixer(self):
        """Test the NHScaleFixer."""

        # Load the molecule
        mol = Molecule('2LCK')

        # Load the sign-inverted data
        path = os.path.dirname(os.path.abspath(__file__))
        data_ref = read_pa_file(path + '/data/2LCK.mr.gz')

        # First, the sign of the N-H RDCs must be fixed. Setup the fixer.
        fixer = SignFixer(mol)

        # The data does not need to be fixed.
        data_fixed, fixes = fixer.fix(data_ref)

        # Setup the fixer
        fixer = NHScaleFixer(mol)

        # The data needs to have the CA-C and C-N+1 RDCs scaled
        data_fixed, fixes = fixer.fix(data_fixed)
        self.assertEqual(len(fixes), 2)

        # Check to see if the overall Q-factor is good
        RMS, Q, data_pred = fixer.fit(data_fixed)

        self.assertLessEqual(Q, 15.0)

    def test_outlierfixer(self):
        """Test the OutlierFixer"""
        # Load the molecule
        mol = Molecule('2MJB')

        # Load the sign-inverted data
        path = os.path.dirname(os.path.abspath(__file__))
        data_ref = read_pa_file(path + '/data/ubq_bicelle_hn-c.pa')

        # Setup the fixer
        fixer = OutlierFixer(mol)

        # The data needs to 10C and 70C removed
        data_fixed, fixes = fixer.fix(data_ref)
        self.assertEqual(len(fixes), 1)

        # The new dataset has 2 fewer points
        self.assertEqual(len(data_ref) - 2,
                         len(data_fixed))

        self.assertIn('A.10C', data_ref)
        self.assertIn('A.70C', data_ref)
        self.assertNotIn('A.10C', data_fixed)
        self.assertNotIn('A.70C', data_fixed)


