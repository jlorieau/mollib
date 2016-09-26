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
    """Collect statistics on Ramachandran angles.
    """

    data_path = mollib.core.settings.ramachandran_dataset_path

    def process_measurement(self, molecule):
        """Process the molecule and return the Ramachandran angle statistics.

        Returns
        -------
        measurement_dict: dict
            The dict produced by :meth:`process_measurement`.

            - key: molecule identifier (str)
            - value: list of tuples of phi-psi angles [(float, float),]
        """
        molecule = super(RamachandranStatistics, self).process_measurement(molecule)

        # Hydrogenate the molecule
        mollib.hydrogens.add_hydrogens(molecule)

        # Get the hydrogen bonds for the molecule
        hbonds = mollib.hbonds.find_hbond_partners(molecule)

        # Group the hbonds into (chain.id, residue.number) and
        #  major_classifications
        residue_classes = dict()

        for hbond in hbonds:
            try:
                donor_residue = hbond.donor.atom2.residue
                acceptor_residue = hbond.acceptor.atom2.residue
            except AttributeError:
                continue
            major_classification = hbond.major_classification
            minor_classification = hbond.minor_classification

            if major_classification == mollib.hbonds.settings.major_bb_bb_amide:
                key = (donor_residue.chain.id, donor_residue.number)
                residue_classes[key] = minor_classification

                key = (acceptor_residue.chain.id, acceptor_residue.number)
                residue_classes[key] = minor_classification

        # Get the Ramachandran angles and group them.
        return_dict = dict()

        for residue in molecule.residues:
            phi_psi = residue.ramachandran_angles
            key = (residue.chain.id, residue.number)
            group = residue_classes.get(key, 'No hydrogen bonds')

            # None values are not saved
            if all(i is not None for i in phi_psi):
                return_dict.setdefault(group, list()).append(phi_psi)

        return return_dict

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

            - key: molecule identifier (str)
            - value: list of tuples of phi-psi angles [(float, float),]
        """
        super(RamachandranStatistics, self).process_data(measurement_dict)

        path = os.path.join(self.root_path, '..', self.data_path)

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
