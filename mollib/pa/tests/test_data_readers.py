"""
Tests for the data reader functions.
"""

import unittest
import os

from mollib.pa.data_readers import read_pa_file, re_mr

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

    def test_mr_regex(self):
        """Test the regex used to match RDC/RACS data in MR format files."""

        # The following comes from 2mjb.mr
        format1 = """
        assign ( resid   500 and name   OO)
               ( resid   500 and name    Z)
               ( resid   500 and name    X)
               ( resid   500 and name    Y)
               ( resid     2 and name    N)
               ( resid     2 and name   HN) -8.1700  0.0000  0.0000
        """

        # Check that the string matched and was properly parsed
        match = re_mr.search(format1)
        self.assertIsNotNone(match)

        d = match.groupdict()
        self.assertEqual(d['coord_num'], '500')
        self.assertEqual(d['res_i'], '2')
        self.assertEqual(d['name_i'], 'N')
        self.assertEqual(d['res_j'], '2')
        self.assertEqual(d['name_j'], 'HN')
        self.assertEqual(d['value_j'], '-8.1700')

        # The following comes from 2oed.mr
        format2 = """
        assign (               resi 501 and name  OO  )
               (               resi 501 and name  Z   )
               (               resi 501 and name  X   )
               (               resi 501 and name  Y   )
               ( segi GB3N and resi   4 and name  C   )
               ( segi GB3N and resi   4 and name  CA  )  -0.906 0.20 0.20
        """

        # Check that the string matched and was properly parsed
        match = re_mr.search(format2)
        self.assertIsNotNone(match)

        d = match.groupdict(); print(d)
        self.assertEqual(d['coord_num'], '501')
        self.assertEqual(d['res_i'], '4')
        self.assertEqual(d['name_i'], 'C')
        self.assertEqual(d['res_j'], '4')
        self.assertEqual(d['name_j'], 'CA')
        self.assertEqual(d['value_j'], '-0.906')

