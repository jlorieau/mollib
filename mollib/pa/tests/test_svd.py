import unittest
import os
import glob

from mollib import Molecule
from mollib.pa import (Process, read_pa_file, calc_pa_SVD, report_tables,
                       find_outliers)
from mollib.pa import settings
from mollib.utils import FormattedStr, MDTable
from mollib.utils.interactions import sort_func


class TestSVD(unittest.TestCase):

    def test_svd_ubiquitin_NH(self):
        """Test the SVD of ubiquitin NH RDCs using PNA data"""

        mol = Molecule('2MJB')
        process = Process(mol)
        magnetic_interactions = process.process()

        # Load the data
        path = os.path.dirname(os.path.abspath(__file__))
        data = read_pa_file(path + '/data/ubq_pna.pa')

        # Fit the data with a SVD
        data_pred, Saupe_components, stats = calc_pa_SVD(magnetic_interactions,
                                                         data)

        # Make sure all of the RDCs are within 6.0Hz
        for label in sorted(data_pred, key=sort_func):
            if label in data:
                self.assertLess(abs(data[label].value - data_pred[label].value),
                                6.0)

    def test_svd_reverse(self):
        """Tests the SVD of dipoles with atom names listed forward and in
        reverse order."""

        # Load the data
        path = os.path.dirname(os.path.abspath(__file__))
        data = read_pa_file(path + '/data/2kxa_sag.inp')

        # Load the reverse data
        data_rev = read_pa_file(path + '/data/2kxa_sag_rev.inp')

        # Process the A-matrix of the molecule
        mol = Molecule('2MJB')
        process = Process(mol)
        process_rev = Process(mol)

        magnetic_interactions = process.process(labels=data.keys())
        magnetic_interactions_rev = process_rev.process(labels=data_rev.keys())

        # Fit the data with a SVD
        data_pred, _, stats = calc_pa_SVD(magnetic_interactions, data)
        data_pred_rev, _, stats_rev = calc_pa_SVD(magnetic_interactions_rev,
                                                  data_rev)

        # Make sure both match
        self.assertEqual(stats['Overall']['Q (%)'],
                         stats_rev['Overall']['Q (%)'])
        self.assertEqual(stats['Overall']['count'],
                         stats_rev['Overall']['count'])

        # Make sure that only 'N-H' is in the stats dict, and not the reverse,
        # 'H-N
        for d in (stats, stats_rev):
            print(d)
            self.assertIn('N-H', d)
            self.assertNotIn('H-N', d)

    def test_svd_ubiquitin_RDC_RACS(self):
        """Test the SVD of ubiquitin with the H-N RDCs and C, N and H RACS from
        Cornilescu JACS 2000."""

        # Load the data
        path = os.path.dirname(os.path.abspath(__file__))
        data = {}
        for filename in glob.glob(path + '/data/ubq_bicelle_*.pa'):

            returned_data = read_pa_file(filename)
            data.update(returned_data)

        # Process the A-matrix of the molecule
        mol = Molecule('2MJB')
        process = Process(mol)
        magnetic_interactions = process.process(labels=data.keys())

        # Fit the data with a SVD
        data_pred, _, stats = calc_pa_SVD(magnetic_interactions, data)

        # Check that the fit stats are reasonable
        self.assertEqual(stats['Overall']['count'], 252)
        self.assertLessEqual(stats['Overall']['Q (%)'], 30.)

        self.assertEqual(stats['N-H']['count'], 63)
        self.assertLessEqual(stats['N-H']['Q (%)'], 15.)
        self.assertLessEqual(stats['N-H']['RMS'], 2.0)

        self.assertEqual(stats['C']['count'], 63)
        self.assertLessEqual(stats['C']['Q (%)'], 18.)
        self.assertLessEqual(stats['C']['RMS'], 9.8)

        self.assertEqual(stats['N']['count'], 63)
        self.assertLessEqual(stats['N']['Q (%)'], 17.)
        self.assertLessEqual(stats['N']['RMS'], 11.5)

        self.assertEqual(stats['H']['count'], 63)
        self.assertLessEqual(stats['H']['Q (%)'], 45.)
        self.assertLessEqual(stats['H']['RMS'], 1.6)

    def test_stats(self):
        """Test the calculation of the Q-factor using PNA data."""
        mol = Molecule('2MJB')

        # Calculate the Q-factor first from the bond lengths
        process = Process(mol)
        settings.calculate_from_bonds = True
        magnetic_interactions = process.process()

        # Load the data
        path = os.path.dirname(os.path.abspath(__file__))
        data = read_pa_file(path + '/data/ubq_pna.pa')

        # Fit the data with a SVD
        data_pred, Saupe_components, stats1 = calc_pa_SVD(magnetic_interactions,
                                                          data)

        # The fit Q-factor should be better than 10%
        self.assertLessEqual(stats1['Overall']['Q (%)'], 10.0)

        # Now calculate it from static DCCs. This Q-factor should be different
        # (yet still low)
        process = Process(mol)
        settings.calculate_from_bonds = False
        magnetic_interactions = process.process()

        # Fit the data with a SVD
        data_pred, Saupe_components, stats2 = calc_pa_SVD(magnetic_interactions,
                                                          data)

        # The fit Q-factor should be better than 10%
        self.assertLessEqual(stats2['Overall']['Q (%)'], 10.0)

        # The RMS should be different between calculating the
        # RDCs using bond lengths vs static values
        self.assertNotEqual(stats1['Overall']['RMS'],
                            stats2['Overall']['RMS'])

    def test_outliers(self):
        """Test the find_outliers function."""
        # Load the molecule and process its interaction arrays.
        mol = Molecule('2MJB')
        process = Process(mol)
        magnetic_interactions = process.process()

        # Load the data
        path = os.path.dirname(os.path.abspath(__file__))
        data = read_pa_file(path + '/data/ubq_pna.pa')

        data_pred, Saupe_components, stats = calc_pa_SVD(magnetic_interactions,
                                                         data)

        warning, bad = find_outliers(data, data_pred)

        self.assertEqual(warning, ['A.11N-H'])
        self.assertEqual(bad, [])

        # Load a second dataset with intentional outliers
        path = os.path.dirname(os.path.abspath(__file__))
        data = read_pa_file(path + '/data/ubq_bicelle_outliers.pa')

        data_pred, Saupe_components, stats = calc_pa_SVD(magnetic_interactions,
                                                         data)

        warning, bad = find_outliers(data, data_pred)

        self.assertEqual(warning, [])
        self.assertEqual(bad, ['A.16N-H', 'A.32N-H'])

    def test_report_tables(self):
        """Test the formatting of the report tables for the SVD fits."""
        # Disable string formatting to more easily match strings
        FormattedStr.color_term = False

        # Load the molecule and process its interaction arrays.
        mol = Molecule('2MJB')
        process = Process(mol)
        magnetic_interactions = process.process()

        # Load the data
        path = os.path.dirname(os.path.abspath(__file__))
        data = read_pa_file(path + '/data/ubq_pna.pa')

        # Fit the data with a SVD
        data_pred, Saupe_components, stats = calc_pa_SVD(magnetic_interactions,
                                                         data)

        # Generate the report tables
        tables = report_tables(data, data_pred)

        # Make sure the fit and predicted tables are present
        for key in ('fit', 'C-N_pred', 'CA-HA_pred', 'N-H_pred', 'C-CA_pred'):
            self.assertIn(key, tables)

            table = tables[key]
            self.assertIsInstance(table, MDTable)
