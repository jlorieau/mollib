"""
Statistics for hydrogen bonds.
"""
import os

import numpy as np

from .statistics import Statistics
from . import settings
import mollib.core
from mollib.hydrogens import add_hydrogens
from mollib.hbonds import find_hbond_partners


class HbondStatistics(Statistics):
    """Collect statistics on hydrogen bonds. A subclass of the
    :class:`mollib.statistics.statistics.Statistics` class.
    """

    def __init__(self, *args, **kwargs):
        super(HbondStatistics, self).__init__(*args, **kwargs)
        self.data_path = mollib.core.settings.hbond_dataset_path

    def process_measurement(self, molecule):
        """Process the molecule and return the Hbond statistics.

        Parameters
        ----------
        molecule: :obj:`mollib.Molecule`
            The molecule to process.

        Returns
        -------
        measurement_dict: dict
            The measurement dict.

            - **key**: classification (type_classification, 
              major_classification, minor_classification), tuple of str
            - **value**: list of dicts
            
              - 'distances': {'a1d1': float, 'a1d2': float, 'a2d1': float,
                'a2d2': float}
              - 'angles': {'theta': float, 'phi': float}
        """
        molecule = super(HbondStatistics, self).process_measurement(molecule)

        # Get the hydrogen bonds for the molecule.
        add_hydrogens(molecule)
        hbonds = find_hbond_partners(molecule)

        # Prepare the measure_dict to return the values
        measure_dict = {}

        for hbond in hbonds:
            # Create the dict key from the major_classification,
            # major_classification and minor_classification. These have to be
            # put together into a string so that json can dumps the dict.
            key = (hbond.type_classification,
                   hbond.major_classification,
                   hbond.minor_classification)
            key = '__'.join(key)
            values_list = measure_dict.setdefault(key, list())

            values_list.append({'distances': hbond.distances,
                                'angles': hbond.angles})

        return measure_dict

    def process_data(self, measurement_dict):
        """Process the output datasets from the data_dict.

        - Creates 1d histogram files for the energies, E(kT), as a function of
          the hbond distances. These are saved as '.npz' files with the
          'd1a1' and 'hist1d' arrays.


        Parameters
        ----------
        measurement_dict: dict
            - **key**: molecule identifier, str
            - **value**: data dict. see output from 
              :meth:`process_measurements`
        """
        super(HbondStatistics, self).process_data(measurement_dict)

        path = self.data_path

        hist_range = (1.0, 3.5)
        bins = 25

        # Convert the measurement_dict into a dict organized by secondary
        # structure classification
        class_dict = {}
        for identifier, return_dict in measurement_dict.items():
            for classification, hbond_dict_list in return_dict.items():
                major, minor, modifier = classification.split('__')
                key = (major, minor)
                d1a1_theta_phi_list = class_dict.setdefault(key, list())
                values = [(d['distances']['d1a1'],
                           d['angles']['theta'],
                           d['angles']['phi']) for d in hbond_dict_list]
                d1a1_theta_phi_list += values

        # Save the 2d histograms numpy arrays
        for classification, d1a1_theta_phi in class_dict.items():
            #d1a1, phi, psi = zip(*d1a1_theta_phi)

            #d1a1 = np.array(d1a1)
            #phi = np.array(phi)
            #psi = np.array(psi)
            d1a1_theta_phi = np.array(d1a1_theta_phi)

            bins = settings.hbond_histogram_bins
            hist, edges = np.histogramdd(d1a1_theta_phi, bins=bins)

            # Convert the histograms into energies (units of kT)
            hist = -1. * np.log(hist + 0.1)
            minimum = np.min(hist)
            hist -= minimum

            # Turn classification tuples into strings separated by '__'
            # characters
            classification = [c for c in classification if c]
            classification = '__'.join(classification)

            filename = os.path.join(path, classification + '.npz')
            np.savez(filename, edges, hist)
