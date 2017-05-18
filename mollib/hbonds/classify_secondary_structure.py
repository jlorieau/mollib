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


# TODO: This file should be split into another plugin, rama


def _group_runs(l, tolerance=2):
    """Generator to group contiguous chain ids and numbers in a list (l) with 
    a tolerance for  skipped numbers.
    
    Adapted from stackoverflow #21142231.
    
    Examples
    --------
    >>> g = _group_runs([('A', 2), ('A', 3), ('A', 4), 
    ...                  ('A', 8), ('A', 10), ('A', 12 )])
    >>> list(g)
    [[('A', 2), ('A', 3), ('A', 4)], [('A', 8), ('A', 10), ('A', 12)]]
    >>> g = _group_runs([('A', 2), ('A', 3), ('A', 4), 
    ...                  ('B', 4), ('B', 5), ('B', 6 )])
    >>> list(g)
    [[('A', 2), ('A', 3), ('A', 4)], [('B', 4), ('B', 5), ('B', 6)]]
    >>> list(_group_runs([]))
    []
    """
    # Sorting is needed to properly group items
    l = sorted(l)

    out = []

    if l:
        last_chain = l[0][0]
        last_item = l[0][1]
        for i, j in l:
            if j - last_item > tolerance or i != last_chain:
                yield out
                out = []
            out.append((i, j))
            last_item = j
            last_chain = i
        yield out


def _is_sheet(residue):
    """Return true if the residue's Ramachandran angles fall within sheet 
    values."""
    phi, psi = residue.ramachandran_angles
    if phi is None or psi is None:
        return False

    # Determine if the phi/psi angles are within a sheet range.
    phi_min, phi_max = settings.beta_phi
    psi_min, psi_max = settings.beta_psi

    return phi_min <= phi <= phi_max and psi_min <= psi <= psi_max


def classify_residues(molecule):
    """Classify the residues of a molecule based on their backbone-
    backbone amide hydrogen bonds.

    Parameters
    ----------
    molecule: :obj:`mollib.Molecule`
        The molecule whose residues are to be classified. The residues gain
        the following attributes:

        - `hbond_classification` with a 'str' to the hbond classification.
        - `hbond_modifier` with a 'str' to the hbond minor classification
           modifier.

    skip_energy: bool, optional
        If True, the ramachandran energy will not be calculated.
    """
    # Hydrogenate the molecule and detect its hydrogen bonds. These are
    # already classified by type
    add_hydrogens(molecule)
    hbonds = find_hbond_partners(molecule)

    # Keep track of groups of classifications, like turns
    turns = {settings.minor_beta_turnI, settings.minor_beta_turnII,
             settings.minor_beta_turnIp, settings.minor_beta_turnIIp}

    # Assign the residue secondary structure based on the hbond
    # classifications
    classification = {}  # {(chain.id, residue.number): classification}
    for hbond in hbonds:
        # Only classify based on backbone-backbone amide hydrogen bonds
        if hbond.major_classification != settings.major_bb_bb_amide:
            continue

        # Classify based on the residue identity for the atom2 (heavy atom)
        # of the acceptor and donor dipoles.
        try:
            donor_res = hbond.donor.atom2.residue
            donor_chain = donor_res.chain.id
            acceptor_res = hbond.acceptor.atom2.residue
            acceptor_chain = acceptor_res.chain.id
        except AttributeError:
            continue

        # Get the hydrogen bond type. The hbond minor classification represents
        # the secondary structure type, like 'sheet'.
        minor_class = hbond.minor_classification
        minor_modifier = ''

        # Isolated hydrogen bonds contain no secondary structure assignment
        # information
        if minor_class in settings.minor_isolated:
            continue

        # Turns are treated specially. If the hbond is for a turn, then in
        # addition to the i and i+3 residues being labeled as turn, so should
        # the i+1 and i+2
        if minor_class in turns:
            res_i = min(donor_res.number, acceptor_res.number)
            res_j = max(donor_res.number, acceptor_res.number)

            # Get the chain id for the donor and acceptor residues
            if donor_chain != acceptor_chain:
                continue
            chain = donor_chain

            # Assign the residues i+1 and i+2 in the turn as turn residues.
            residue_nums = range(res_i + 1, res_j)
            for i in residue_nums:
                classification[(chain, i)] = (minor_class, minor_modifier)
            continue

        # Assign both the donor and acceptor residues
        for count, res in enumerate((donor_res, acceptor_res)):
            # Get the chain id for the donor and acceptor residues
            if donor_chain != acceptor_chain:
                continue
            chain = donor_chain

            # Some minor modifiers do not pertain to the donor
            # residue, where count==0, like the N-terminal residues of
            # helices
            if count == 1 and hbond.minor_modifier == settings.minor_N:
                minor_modifier = hbond.minor_modifier

            # Some minor modifiers do not pertain to the acceptor
            # residue, where count==0, like the C-terminal residues of
            # helices
            if count == 0 and hbond.minor_modifier == settings.minor_C:
                minor_modifier = hbond.minor_modifier

            # If the residue is already assigned and it isn't 'isolated'
            # then skip it
            if (res.number in classification and
                classification[(chain, res.number)][0] !=
               settings.minor_isolated):
                continue

            classification[(chain, res.number)] = (minor_class, minor_modifier)

    # Fill in the classifications. Some classification, like sheets, may pertain
    # to all residues within a contiguous block of sheet residues.
    sheet_types = {settings.minor_beta}
    sheet_residues = [k for k, v in classification.items()
                      if v[0] in sheet_types]

    # Find the contiguous sheet residue numbers and assign sheet residues
    for group in _group_runs(sheet_residues, tolerance=2):
        if len(group) == 0 :
            continue

        # Determine the first and last residue in the sheet group
        first_chain, first_res_num = group[0]
        last_chain, last_res_num = group[-1]

        # Skip if the sheet consists of residues from different chains
        if first_chain != last_chain:
            continue

        # Go through the residues in the group, and make sure they are sheets
        # Start with the residue before the last_chain, and the residue after
        # the first chain.
        for i in range(first_res_num - 1, last_res_num + 2):
            try:
                residue = molecule[first_chain][i]
            except KeyError:
                continue

            # Determine whether the residue counts as a beta sheet
            if not _is_sheet(residue):
                continue

            # Overwrite existing classifications, if it hasn't be assigned yet.
            key = (first_chain, i)
            if key not in classification:
                classification[key] = (settings.minor_beta, '')

    # Classify the residues. The classification dict has been populated
    # {residue.number(int): classification(str}
    # Now use it to classify the residues
    for residue in molecule.residues:
        try:
            chain = residue.chain.id
            res_num = residue.number
        except AttributeError:
            continue

        res_class = classification.get((chain, res_num), ('', ''))

        residue.classification = res_class

        # Add the ramachandran energy attribute to the residue
        add_energy_ramachandran(residue)


#: The Ramachandran energy datasets
energy_ramachandran_datasets = None


def add_energy_ramachandran(residue):
    """Add the Ramachandran energy (kT) to the residue.

    The energy is assigned based on the `hbond_classification` and the energies
    computed from the statistics module data of RamachandranStatistics.

    Parameters
    ----------
    residue: :obj:`mollib.Residue`
        The residue gains the following attributes:

    Added Residues Attributes
    -------------------------
    energy_ramachandran: float
        The Ramachandran energy (in kT)
    ramachandran_dataset: str
        The name of the dataset used for the Ramachandran energy.
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

    # Get the energy for this residue's classification
    classification = residue.classification
    if (classification is None or
        not classification[0] or
       classification[0] == mollib.hbonds.settings.minor_isolated):

        # The 'No Hydrogen Bonds' and 'isolated' classification is the same as
        # the 'No classification'
        res_class = 'No classification'
        modifier = ''
    else:
        res_class = (classification[0] if classification is not None
                     else 'No classification')
        modifier = classification[1] if classification is not None else ''

    # Some energies are classified by hbond_classification and hbond_modifier.
    # These have the name of both, separated by '__'
    full_res_class = '__'.join((res_class, modifier))

    if residue.name == 'GLY' and 'Gly' in energy_ramachandran_datasets:
        # Glycines are treated specially because they are more flexible in
        # their backbone torsion angles
        phi_1d, psi_1d, energy_2d = energy_ramachandran_datasets['Gly']
        residue.ramachandran_dataset = 'Gly'

    elif full_res_class in energy_ramachandran_datasets:
        # First, try to find the corresponding dataset using the
        # hbond_classification and hbond_modifier
        phi_1d, psi_1d, energy_2d = energy_ramachandran_datasets[full_res_class]
        residue.ramachandran_dataset = full_res_class

    elif res_class in energy_ramachandran_datasets:
        # Next, try to find the corresponding dataset using just the
        # hbond_classification
        phi_1d, psi_1d, energy_2d = energy_ramachandran_datasets[res_class]
        residue.ramachandran_dataset = res_class
    elif 'No classification' in energy_ramachandran_datasets:
        (phi_1d, psi_1d,
         energy_2d) = energy_ramachandran_datasets['No classification']
        residue.ramachandran_dataset = 'No classification'
    else:
        return None

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
    except (KeyError, IndexError) as e:
        msg = ("The phi/psi angles ({},{}) for not be found in the "
               "Ramachandran dataset.")
        print(msg.format(phi, psi))
        raise e

    # Assign the energy
    residue.energy_ramachandran = energy
