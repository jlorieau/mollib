"""
Classify residues based on hydrogen bonds.
"""
import glob
import os

import numpy as np

import mollib.core.settings
from mollib.hydrogens import add_hydrogens
from .hbonds import find_hbond_partners
from . import settings


def classify_residues(molecule):
    """Classify the residues of a molecule based on their backbone-
    backbone amide hydrogen bonds.

    Parameters
    ----------
    molecule: :obj:`mollib.Molecule`
        The molecule whose residues are to be classified. The residues gain
        the following attributes:

        - `hbond_classification` with a 'str' to the hbond classification.
    """
    # TODO: make the classification a residue parameter (like molecule
    # parameters) using a parameter mixin class

    # Hydrogenate the molecule and detect its hydrogen bonds. These are
    # already classified by type
    add_hydrogens(molecule)
    hbonds = find_hbond_partners(molecule)

    # Assign the residue secondary structure based on the hbond
    # classifications
    classification = {}  # {residue.number: classification}
    for hbond in hbonds:
        # Only classify based on backbone-backbone amide hydrogen bonds
        if hbond.major_classification != mollib.hbonds.settings.major_bb_bb_amide:
            continue

        # Classify based on the residue identity for the atom2 (heavy atom)
        # of the acceptor and donor dipoles.
        try:
            donor_res = hbond.donor.atom2.residue
            acceptor_res = hbond.acceptor.atom2.residue
        except AttributeError:
            continue

        # Assign both the donor and acceptor residues
        for res in (donor_res, acceptor_res):
            # If not assigned, assigned the classification
            if res.number not in classification:
                minor_class = hbond.minor_classification
                classification[res.number] = minor_class

            # If already assigned, only replace if the current classification
            # is for an 'isolated' hydrogen bond
            else:
                current_class = classification[res.number]
                if current_class == settings.minor_isolated:
                    minor_class = hbond.minor_classification
                    classification[res.number] = minor_class

    # Classify the residues. The classification dict has been populated
    # {residue.number(int): classification(str}
    # Now use it to classify the residues
    for residue in molecule.residues:
        res_class = classification.get(residue.number, '')
        residue.hbond_classification = res_class

        # Add the ramachandran energy attribute to the residue
        add_ramachandran_energy(residue)


#: The Ramachandran energy datasets
energy_ramachandran_datasets = None


def add_ramachandran_energy(residue):
    """Add the Ramachandran energy (kT) to the residue.

    The energy is assigned based on the `hbond_classification` and the energies
    computed from the statistics module data of RamachandranStatistics.

    Parameters
    ----------
    residue: :obj:`mollib.Residue`
        The residue gains the following attributes:

        - `energy_ramachandran` with a  float to the Ramachandran energy.
    """
    global energy_ramachandran_datasets

    # Load the energy datasets.
    if energy_ramachandran_datasets is None:
        energy_ramachandran_datasets = {}

        path = mollib.core.settings.ramachandran_dataset_path
        path = os.path.join(path, '*.npz')
        filepaths = glob.iglob(path)

        for filepath in filepaths:
            path, filename = os.path.split(filepath)
            name, ext = os.path.splitext(filename)

            # Load the dataset from a numpyz file
            # The datasets are histogram2d datasets with a phi_1d array of
            # the phi angles, a psi_1d array of the psi angles and a
            # energy_2d matrix with the corresponding 2d energy at the
            # phi/psi angle.
            try:
                arrays = np.load(filepath)
                phi_1d, psi_1d, energy_2d = (arrays['arr_0'],
                                             arrays['arr_1'],
                                             arrays['arr_2'])
            except:
                continue

            energy_ramachandran_datasets[name] = (phi_1d, psi_1d, energy_2d)

        # The 'No Hydrogen Bonds' classification is the same as the ''
        # classification above
        if 'No hydrogen bonds' in energy_ramachandran_datasets:
            v = energy_ramachandran_datasets['No hydrogen bonds']
            energy_ramachandran_datasets[''] = v

    # Get the energy for this residue's classification
    res_class = getattr(residue, 'hbond_classification', '')

    if res_class in energy_ramachandran_datasets:
        phi_1d, psi_1d, energy_2d = energy_ramachandran_datasets[res_class]
    else:
        phi_1d, psi_1d, energy_2d = energy_ramachandran_datasets['']

    # Find the appropriate phi and psi angles for the residue.
    phi, psi = residue.ramachandran_angles

    # Assigning the energy is only possible if both phi and psi angles are
    # specified
    if not isinstance(phi, float) or not isinstance(psi, float):
        residue.energy_ramachandran = 0.0
        return None

    # Find the corresponding index numbers for the 2d array.
    # The phi_1d array holds the phi angles, the y array holds the psi angle
    # and the z array is a 2d with the energy for the phi/psi angle
    try:
        phi_index = [i for i,p in enumerate(phi_1d) if p >= phi]
        phi_index = phi_index[0]

        psi_index = [i for i, p in enumerate(psi_1d) if p >= psi]
        psi_index = psi_index[0]

        # Get the energy from the energy_2d
        energy = energy_2d[psi_index - 1][phi_index - 1]
    except KeyError as e:
        msg = ("The phi/psi angles ({},{}) for not be found in the "
               "Ramachandran dataset.")
        print(msg.format(phi, psi))
        raise e

    # If avaliable, assign the overall energy too.

    # Assign the energy
    residue.energy_ramachandran = energy







