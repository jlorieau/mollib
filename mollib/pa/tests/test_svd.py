import unittest
import os
import glob

from mollib import Molecule
from mollib.hydrogens import add_hydrogens
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
        mol = Molecule('2KXA')
        process = Process(mol)
        process_rev = Process(mol)

        magnetic_interactions = process.process(labels=data.keys())
        magnetic_interactions_rev = process_rev.process(labels=data_rev.keys())

        # Fit the data with a SVD
        data_pred, _, stats = calc_pa_SVD(magnetic_interactions, data)
        data_pred_rev, _, stats_rev = calc_pa_SVD(magnetic_interactions_rev,
                                                  data_rev)

        # Make sure both match and that the Q-factors are better than 15%
        self.assertEqual(stats['Overall']['Q (%)'],
                         stats_rev['Overall']['Q (%)'])
        self.assertLessEqual(stats['Overall']['Q (%)'], 15.0)
        self.assertLessEqual(stats_rev['Overall']['Q (%)'], 15.0)
        self.assertEqual(stats['Overall']['count'],
                         stats_rev['Overall']['count'])

        # Make sure that only 'N-H' is in the stats dict, and not the reverse,
        # 'H-N
        for d in (stats, stats_rev):
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

    def test_multi_conformer(self):
        """Test the SVD of a multiconformer refinement with 2LWA."""

        # Load the molecules
        path = os.path.dirname(os.path.abspath(__file__))
        mol1 = Molecule(path + '/data/2lwa_struc-a.pdb')
        mol2 = Molecule(path + '/data/2lwa_struc-b.pdb')
        mol3 = Molecule(path + '/data/2lwa_struc-c.pdb')

        # Load the data
        data = read_pa_file(path + '/data/2lwa.inp')

        # Flip the sign of N-H couplings
        for k, v in data.items():
            if k.endswith('N-H'):
                v.value *= -1.

        # Process the A-matrix of the molecule
        process = Process([mol1, mol2, mol3])
        magnetic_interactions = process.process(labels=data.keys())

        # Fit the data with a SVD
        data_pred, _, stats = calc_pa_SVD(magnetic_interactions, data)

        # Check the statistics
        self.assertLessEqual(stats['N-H']['Q (%)'], 15.)
        self.assertLessEqual(stats['CA-HA']['Q (%)'], 20.)


    def test_multi_conformer_same_structure(self):
        """Test the fit of multiple conformers of the same structure.
        
        This test is used to see if the Saupe matrix and Q-factor values are
        comparable between the single structure fit and the multi-structure
        fit. However, the exact same structure cannot be used in the SVD, as
        this produces undefined weights. Therefore, we will use two similar
        ubiquitin structures: 2MJB and 1UBQ. The 1UBQ structure fits 
        (relatively) poorly to the RDC/RACS data, so its Da is much smaller
        than that of 2MJB, which was refined against RDCs/RACSs.
        """
        # Load the molecules. 1UBQ has to be hydrogenated
        mol_2mjb = Molecule('2MJB')
        mol_1ubq = Molecule('1UBQ')
        add_hydrogens(mol_1ubq, strip=True)

        # Load the data
        path = os.path.dirname(os.path.abspath(__file__))
        data = read_pa_file(path + '/data/ubq_bicelle_hn-c.pa')

        # Process the A-matrix for 2MJB and 2MJB+1UBQ
        process_2mjb = Process(mol_2mjb)
        magnetic_interactions_2mjb = process_2mjb.process(labels=data.keys())

        process_both = Process([mol_2mjb, mol_1ubq])
        magnetic_interactions_both = process_both.process(labels=data.keys())

        # Conduct the SVD for each
        _, _, stats_2mjb = calc_pa_SVD(magnetic_interactions_2mjb, data)
        _, _, stats_both = calc_pa_SVD(magnetic_interactions_both, data)

        # Compare the statistics and the Saupe matrices for both. The 2MJB
        # structure is the first structure in the two structure (both) SVD fit.
        # First, both should have reasonable Q-factors
        self.assertLessEqual(stats_2mjb['Overall']['Q (%)'],
                             15.0)
        self.assertLessEqual(stats_both['Overall']['Q (%)'],
                             17.0)

        # In the SVD of both, the Da for 2MJB should be much larger than 1UBQ
        self.assertGreaterEqual(abs(stats_2mjb['N-H']['Da (Hz)']),
                                12.0)
        self.assertGreaterEqual(abs(stats_both['N-H (1)']['Da (1) (Hz)']),
                                12.0)
        self.assertLessEqual(abs(stats_both['N-H (2)']['Da (2) (Hz)']),
                             2.5)

        # The Da should be within 5% for the 2MJB SVD and the 2MJB
        # structure of the SVD with both structures
        Da_1 = stats_2mjb['N-H']['Da (Hz)']
        Da_2 = stats_both['N-H (1)']['Da (1) (Hz)']
        Da_diff = (Da_1 - Da_2) * 2. / (Da_1 + Da_2)
        self.assertLessEqual(abs(Da_diff), 0.05)

        # The Rh should be within 12%
        Rh_1 = stats_2mjb['N-H']['Rh']
        Rh_2 = stats_both['N-H (1)']['Rh']
        Rh_diff = (Rh_1 - Rh_2) * 2. / (Rh_1 + Rh_2)
        self.assertLessEqual(abs(Rh_diff), 0.12)

        # The Saupe matrix should be within 7.5% as well. However, the smallest
        # component has a relatively large error, so these are compared to the
        # largest component (Szz)
        Szz = stats_2mjb['Saupe']['Szz']
        msg = "The '{}' parameters are {} {} ({}% different)"
        for parameter in ('Szz', 'Sxx', 'Syy'):
            value1 = stats_2mjb['Saupe'][parameter]
            value2 = stats_both['Saupe (1)'][parameter]
            value_diff = (value1 - value2) / Szz
            self.assertLessEqual(abs(value_diff), 0.075,
                                 msg.format(parameter, value1, value2,
                                            value_diff * 100.))

        # The angles should be within 10 degrees as well.
        msg = "The '{}' parameters are {} {} ({} different)"
        for parameter in ("Z (deg)", "Y' (deg)" , "Z'' (deg)"):
            value1 = stats_2mjb['Angles'][parameter]
            value2 = stats_both['Angles (1)'][parameter]
            value_diff = (value1 - value2)
            self.assertLessEqual(abs(value_diff), 10.,
                                 msg.format(parameter, value1, value2,
                                            value_diff))

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
        for key in ('fit', 'pred'):
            self.assertIn(key, tables)

            table = tables[key]
            self.assertIsInstance(table, MDTable)
