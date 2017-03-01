import unittest
import os

from mollib import Molecule
from mollib.pa.process_molecule import Process
from mollib.pa.data_readers import read_pa_file
from mollib.pa.svd import calc_pa_SVD
from mollib.pa.utils import sort_key
from mollib.pa import settings


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
        self.assertLessEqual(stats1['Q'], 0.10)

        # Now calculate it from static DCCs. This Q-factor should be different
        # (yet still low)
        process = Process(mol)
        settings.calculate_from_bonds = False
        magnetic_interactions = process.process()

        # Fit the data with a SVD
        data_pred, Saupe_components, stats2 = calc_pa_SVD(magnetic_interactions,
                                                          data)

        # The fit Q-factor should be better than 10%
        self.assertLessEqual(stats2['Q'], 0.10)

        # The residual sum squared should be different between calculating the
        # RDCs using bond lengths vs static values
        self.assertNotEqual(stats1['RSS'], stats2['RSS'])
