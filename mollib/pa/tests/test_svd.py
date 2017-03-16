import unittest
import os

from mollib import Molecule
from mollib.pa import (Process, read_pa_file, calc_pa_SVD, report_tables,
                       find_outliers, sort_key)
from mollib.pa import settings
from mollib.utils import FormattedStr, MDTable


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
        for label in sorted(data_pred, key=sort_key):
            if label in data:
                self.assertLess(abs(data[label].value - data_pred[label].value),
                                6.0)

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
        self.assertLessEqual(stats1['Q-factor (%)'], 10.0)

        # Now calculate it from static DCCs. This Q-factor should be different
        # (yet still low)
        process = Process(mol)
        settings.calculate_from_bonds = False
        magnetic_interactions = process.process()

        # Fit the data with a SVD
        data_pred, Saupe_components, stats2 = calc_pa_SVD(magnetic_interactions,
                                                          data)

        # The fit Q-factor should be better than 10%
        self.assertLessEqual(stats2['Q-factor (%)'], 10.0)

        # The residual sum squared should be different between calculating the
        # RDCs using bond lengths vs static values
        self.assertNotEqual(stats1['RSS'], stats2['RSS'])

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
        for key in ('fit', 'N-C_pred', 'CA-HA_pred', 'N-H_pred', 'CA-C_pred'):
            self.assertIn(key, tables)

            table = tables[key]
            self.assertIsInstance(table, MDTable)
