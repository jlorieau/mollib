"""
Statistics for Ramachandran angles
"""

import os

import numpy as np

import mollib.core
import mollib.hbonds
import mollib.hydrogens
from .statistics import Statistics
from . import settings


class RamachandranStatistics(Statistics):
    """Collect statistics on Ramachandran angles. A subclass of the
    :class:`mollib.statistics.statistics.Statistics` class.
    """

    def __init__(self, *args, **kwargs):
        super(RamachandranStatistics, self).__init__(*args, **kwargs)
        self.data_path = mollib.core.settings.ramachandran_dataset_path

    def process_measurement(self, molecule):
        """Process the molecule and return the Ramachandran angle statistics.

        Parameters
        ----------
        molecule: :obj:`mollib.Molecule`
            The molecule to process.

        Returns
        -------
        measurement_dict: dict
            The dict produced by :meth:`process_measurement`.

            - **key**: residue classification, str
            - **value**: list of tuples of phi-psi angles, [(float, float),]
        """
        molecule = super(RamachandranStatistics,
                         self).process_measurement(molecule)

        # Get the hydrogen bonds for the molecule. This function also adds
        # the residue.hbond_classification and residue.hbond_modifier
        # attributes
        mollib.hbonds.classify_residues(molecule)

        # Prepare the measure_dict to return the values
        measure_dict = {}

        for residue in molecule.residues:
            phi_psi = residue.ramachandran_angles

            if residue.name == 'GLY':
                major_classification = 'Gly'
                minor_classification = ''
            else:
                classification = residue.classification

                # Group unclassified and isolated hydrogen bond residues
                # together under 'No classification'. Otherwise just group them
                # by their classification.
                if (classification is None or
                    not classification[0] or
                    classification[0] == mollib.hbonds.settings.major_isolated):
                    major_classification = 'No classification'
                    minor_classification = ''
                else:
                    major_classification = classification[0]
                    minor_classification = classification[1]

            # None values are not saved
            if all(i is not None for i in phi_psi):
                # Save the classification

                name = (major_classification if not minor_classification else
                        '__'.join((major_classification, minor_classification)))
                measure_dict.setdefault(name, list()).append(phi_psi)

        return measure_dict

    def process_data(self, measurement_dict):
        """Process the output datasets from the data_dict.

        - Creates 2d histogram files for the energies, E(kT), as a function of
          the Ramachandran angles. These are saved as '.npz' files with the
          'phi', 'psi' and 'hist2d' arrays.
        - Creates background contour plot files (pdf) for the Ramachandran
          probability densities. Each contour represents one kT unit.

        Parameters
        ----------
        measurement_dict: dict
            The dict produced by :meth:`process_measurement`.

            - **key**: molecule identifier, str
            - **value**: list of tuples of phi-psi angles, [(float, float),]
        """
        super(RamachandranStatistics, self).process_data(measurement_dict)

        path = self.data_path

        # Convert the measurement_dict into a dict organized by secondary
        # structure classification
        class_dict = {}
        for identifier, return_dict in measurement_dict.items():
            for classification, phi_psi_list in return_dict.items():
                l = class_dict.setdefault(classification, list())
                phi_psi_list = [(i, j) for i, j in phi_psi_list if
                                isinstance(i, float) and isinstance(j, float)]
                l.extend(phi_psi_list)

        # Save the 2d histograms numpy arrays
        for classification, phi_psi in class_dict.items():
            phi, psi = zip(*phi_psi)

            phi = np.array(phi)
            psi = np.array(psi)
            bins = settings.ramachandran_histogram_bins
            hist2d, phi, psi = np.histogram2d(psi, phi, bins=bins,
                                range=np.array([(-180., 180.), (-180., 180.)]))

            # Convert the histograms into energies (units of kT)
            hist2d = -1. * np.log(hist2d + 0.1)
            minimum = np.min(hist2d)
            hist2d -= minimum

            filename = os.path.join(path, classification + '.npz')
            np.savez(filename, phi, psi, hist2d)
