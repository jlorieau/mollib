"""
Tests for the data reader functions.
"""

import unittest
import os

from mollib.pa.data_readers import read_pa_file

# TODO: Test file reading for interactions with atom names that are backwards
class TestDataReader(unittest.TestCase):

    def test_pa_file(self):
        """Test the reading of the standard pa file."""

        # This dataset has 126 RDCs and RACS
        path = os.path.dirname(os.path.abspath(__file__))
        data = read_pa_file(path + '/data/ubq_bicelle_hn-c.pa')
        self.assertEqual(len(data), 126)

        # Check one RDC and one RACS
        self.assertEqual(data['A.2N-H'].value, 13.2)
        self.assertEqual(data['A.1C'].value, 9.1)

    def test_dc_file(self):
        """Test the reading for the DC format in NMRPipe."""

        # This dataset has 42 RDCs and RACS
        path = os.path.dirname(os.path.abspath(__file__))
        data = read_pa_file(path + '/data/2kxa_sag.inp')
        self.assertEqual(len(data), 42)

        # Check one RDC of each interaction type
        self.assertEqual(data['A.3N-H'].value, -10.5)
        self.assertEqual(data['A.1CA-HA#'].value, -8.3)

    def test_mr_file(self):
        """Test the reading for the MR format used in the PDB."""

        # This dataset has 58 RDCs
        path = os.path.dirname(os.path.abspath(__file__))
        data = read_pa_file(path + '/data/2kxa.mr.gz')
        self.assertEqual(len(data), 58)

        # Check one RDC of each interaction type
        self.assertEqual(data['A.3N-H'].value, 10.5)
        self.assertEqual(data['A.18CA-HA'].value, -17.9)
