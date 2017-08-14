"""
Classify residues based on hydrogen bonds. Additionally, there is a fill_gap
function that is used to identify contiguous secondary structure units that may
not have hydrogen bonds with all residues.
"""
import glob
import os

import numpy as np

import mollib.core.settings
from mollib.hydrogens import add_hydrogens
from .hbonds import find_hbond_partners
from . import settings


# TODO: Make 'Gly' classification


def classify_residues(molecule):
    """Classify the residues of a molecule based on their backbone-
    backbone amide hydrogen bonds.

    Parameters
    ----------
    molecule: :obj:`mollib.Molecule`
        The molecule whose residues are to be classified. 
        
        .. note:: The residue objects (:obj:`mollib.Residue`) has the 
            :attr:`classification` tuple assigned.
        
            - 'classification': (`major_classification`, `minor_classification`)

    skip_energy: bool, optional
        If True, the ramachandran energy will not be calculated.
    """
    # Hydrogenate the molecule and detect its hydrogen bonds. These are
    # already classified by type
    add_hydrogens(molecule)
    hbonds = find_hbond_partners(molecule)

    # Keep track of groups of classifications, like turns
    turns = {settings.major_beta_turnI, settings.major_beta_turnII,
             settings.major_beta_turnIp, settings.major_beta_turnIIp}

    # Assign the residue secondary structure based on the hbond
    # classifications
    classification = {}  # {(chain.id, residue.number): classification}
    for hbond in hbonds:
        # Only classify based on backbone-backbone amide hydrogen bonds
        if hbond.type_classification != settings.type_bb_bb_amide:
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

        # Get the hydrogen bond major classification. The hbond major
        # classification represents the secondary structure type, like 'sheet'.
        major_class = hbond.major_classification

        # Isolated hydrogen bonds contain no secondary structure assignment
        # information. Skip these.
        if major_class in settings.major_isolated:
            continue

        # Turns are treated specially. If the hbond is for a turn, residues
        # i+1 and i+2 shouuld be labeled as a turn (even though the i/i+4 is
        # hydrogen bond)
        if major_class in turns:
            res_i = min(donor_res.number, acceptor_res.number)
            res_j = max(donor_res.number, acceptor_res.number)

            # Get the chain id for the donor and acceptor residues
            if donor_chain != acceptor_chain:
                continue
            chain = donor_chain

            # Assign the residues i+1 and i+2 in the turn as turn residues.
            residue_nums = range(res_i + 1, res_j)
            for i in residue_nums:
                classification[(chain, i)] = (major_class, '')
            continue

        # Assign both the donor and acceptor residues. This code block assigns
        # the major_classification of residues. The minor_classification is
        # later labeled by the assign_block function.
        for count, res in enumerate((donor_res, acceptor_res)):
            # Reset the minor classification to nothing
            minor_class = ''

            # Get the chain id for the donor and acceptor residues
            if donor_chain != acceptor_chain:
                continue
            chain = donor_chain

            # # Some minor classifications do not pertain to the donor
            # # residue, where count==0, like the N-terminal residues of
            # # helices
            # if count == 1 and hbond.minor_classification == settings.minor_N:
            #     minor_class = hbond.minor_classification
            #
            # # Some minor classifications do not pertain to the acceptor
            # # residue, where count==0, like the C-terminal residues of
            # # helices
            # if count == 0 and hbond.minor_classification == settings.minor_C:
            #     minor_class = hbond.minor_classification

            # If the residue is already assigned and it isn't 'isolated'
            # then skip it
            if (res.number in classification and
                classification[(chain, res.number)][0] !=
               settings.major_isolated):
                continue

            classification[(chain, res.number)] = (major_class, minor_class)

    # Find the contiguous blocks of secondary structure elements from the
    # hydrogen bonds and fill in gaps in the primary sequence
    if settings.assign_blocks:
        # Determine whether existing assignments should be overwritten
        overwrite_assignments = settings.assign_blocks_overwrite

        # Assign blocks for alpha-helices
        if settings.assign_blocks_alpha:
            #  Get the parameters to fill the gaps
            kwargs = {
                'extend_termini': settings.assign_blocks_alpha_extend_termini,
                'label_N_term': settings.assign_blocks_alpha_label_N_term,
                'label_C_term': settings.assign_blocks_alpha_label_C_term,
                'gap_tolerance': settings.assign_blocks_alpha_gap_tolerance,
                'overwrite_assignments': overwrite_assignments,
                      }

            assign_blocks(molecule, classification, settings.major_alpha,
                          _is_helix, **kwargs)

        # Assign blocks and fill gaps and  for beta-strands and sheets
        if settings.assign_blocks_beta:
            # Get the parameters to fill the gaps
            kwargs = {
                'extend_termini': settings.assign_blocks_beta_extend_termini,
                'label_N_term': settings.assign_blocks_beta_label_N_term,
                'label_C_term': settings.assign_blocks_beta_label_C_term,
                'gap_tolerance': settings.assign_blocks_beta_gap_tolerance,
                'overwrite_assignments': overwrite_assignments,
                      }

            assign_blocks(molecule, classification, settings.major_beta,
                          _is_sheet, **kwargs)

        # Assign blocks and fill gaps for 310-helices
        if settings.assign_blocks_310:
            # 310-helices are typically 4-5 residues long. For a 4-residue
            # 310-helix, the i and i+3 residues are hydrogen bonded and labeled
            # as 310-helix--but the i+1/i+2 are not. This will fill in that gap

            # Get the parameters to fill the gaps
            kwargs = {
                'extend_termini': settings.assign_blocks_310_extend_termini,
                'label_N_term': settings.assign_blocks_310_label_N_term,
                'label_C_term': settings.assign_blocks_310_label_C_term,
                'gap_tolerance': settings.assign_blocks_310_gap_tolerance,
                'overwrite_assignments': overwrite_assignments,
                      }

            assign_blocks(molecule, classification, settings.major_310,
                          _is_helix, **kwargs)

        # TODO: add gap filling for pi-helices and potentiall alpha-helices.

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


# The cached Ramachandran energy datasets
energy_ramachandran_datasets = None


def add_energy_ramachandran(residue):
    """Add the Ramachandran energy (kT) to the residue.

    The energy is assigned based on the `hbond_classification` and the energies
    computed from the statistics module data of RamachandranStatistics.

    Parameters
    ----------
    residue: :obj:`mollib.Residue`
        The residue gains the following attributes

        .. note:: The following attributes are added to the residue:
        
            - 'energy_ramachandran': The Ramachandran energy (in kT), float  
            - 'ramachandran_dataset': The name of the dataset used for the 
              Ramachandran energy, str
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
       classification[0] == mollib.hbonds.settings.major_isolated):

        # The 'No Hydrogen Bonds' and 'isolated' classification is the same as
        # the 'No classification'
        res_class = 'No classification'
        minor_class = ''
    else:
        res_class = (classification[0] if classification is not None
                     else 'No classification')
        minor_class = classification[1] if classification is not None else ''

    # Some energies are classified by hbond_classification and hbond_modifier.
    # These have the name of both, separated by '__'
    full_res_class = '__'.join((res_class, minor_class))

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
        residue.ramachandran_dataset = 'No classification'
        residue.energy_ramachandran = 0.0
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


def assign_blocks(molecule, classifications, classification_type, dihedral_test,
                  extend_termini=False, label_N_term=0, label_C_term=0,
                  gap_tolerance=1, overwrite_assignments=False):
    """Assign blocks of secondary structure units, including those with gaps.
    A secondary structure block consists of a contiguous series of residues with
    the same secondary structure major_classification.
    
    This function also assigns the minor 'N-term' and 'C-term' classifications.
    
    Gaps occur in the secondary structure assignment from hydrogen bonds, for 
    example, with beta-strands on the edges of beta sheets. This function finds
    stretches of secondary structure assignments, it checks the dihedral angles
    and fills in gaps in the stretches. For a sheet: 'E E E E' becomes 
    'EEEEEEE'.
    
    Parameters
    ----------
    molecule: :obj:`mollib.Molecule`
        The molecule object to classify the secondary structure elements.
    classifications: dict
        A dict with the classifications.
        
          - **key**: (chain.id, residue.number). ex: ('A', 31)
          - **value**: (major_classification, major_classification).
            ex: ('alpha-helix', 'N-term')
    classification_type: str
        The name of the major classification to check. ex: 'alpha-helix'
    dihedral_test: function or None
        - A test function that takes a :obj:`mollib.Residue` and returns True
          if the residue's dihedral angles are within range for the 
          'classification_type'.
        - If None is specified, then the dihedral angles of residues will not
          be tested.
    extend_termini: bool, optional
        If True, the previous and subsequence residues of each contiguous
        stretch of residue classification will be checked to see if they fall
        within the dihedral angle range as well.
    label_N_term: int, optional
        Label the first 'N' number of residues in the contiguous block as 
        'N-term'
    label_C_term: int, optional
        Label the last 'N' number of residues in the contiguous block as 
        'C-term'
    gap_tolerance: int, optional
        The assignment of contiguous stretches of a secondary structure 
        assignment will tolerate this number of 'gaps' in the residue numbers.
        For a gap_toleranace of 2 and a checked sheet assignment, the following
        group 'E E E E' will be treated as a single contiguous block of sheet
        assignments. gap_tolerance should be greater than 0.
    overwrite_assignments: bool, optional
        If True, classification assignments will be overwritten, if an
        assignments has already been made for a given residue.

    Returns
    -------
    None
    """
    assert gap_tolerance > 0

    # Filter the classifications (keys) that match the classification_type
    filtered_keys = [k for k, v in classifications.items()
                     if v[0] == classification_type]

    # Group the keys by contiguous blocks
    groups = _group_runs(filtered_keys, tolerance=gap_tolerance)

    # Go through the groups and make the secondary structure classification
    # assignments
    for group in groups:
        # If the group has no items, skip it
        if len(group) == 0:
            continue

        # Determine the chain id and first and last residue of the group
        first_chain, first_res_num = group[0]
        last_chain, last_res_num = group[-1]

        # If the residues do not belong to the same chain, then there is a
        # bug in the assignment of the contiguous block
        assert first_chain == last_chain

        # Determine the contiguous residue numbers (for the given chain)
        # to evaluate whether they are part of the classification_type
        first = min(first_res_num, last_res_num)
        last = max(first_res_num, last_res_num)

        if extend_termini:
            # If extend terminii, check the previous and subsequence residues
            # as well as part of the secondary structure element
            residue_numbers = range(first - 1, last + 2)
        else:
            # Otherwise, check all residues between the first and last residues
            # (inclusive)
            residue_numbers = range(first, last + 1)

        # Go through the residue numbers and see if they should be added to the
        # contiguous stretch of secondary structure assignments. We first
        # make a contiguous_keys list for the residues that belong to the
        # contiguous stretch. These will subsequently be assigned to the
        # classification

        contiguous_keys = []
        for count, i in enumerate(residue_numbers):
            # If needed, check to see if the residue falls within the dihedral
            # angle ranges with the dihedral_test function.
            if dihedral_test is not None:
                try:
                    residue = molecule[first_chain][i]
                except KeyError:
                    continue
                if not dihedral_test(residue):
                    continue

            # At this point, the residue can be assigned.
            key = (first_chain, i)
            contiguous_keys.append(key)

        # Now go through the contiguous_keys and assign them to the
        # classifications. This has to be done separately so that the 'N-term'
        # and 'C-term' modifiers are properly assigned
        major_classification = classification_type
        for count, key in enumerate(contiguous_keys):
            current_classification = classifications.get(key, None)

            # Determine whether the 'N-term' or 'C-term' should be labeled
            if count < label_N_term:
                minor_classification = 'N-term'
            elif count >= len(contiguous_keys) - label_C_term:
                minor_classification = 'C-term'
            else:
                minor_classification = ''

            # If the current_classification is not assigned (i.e. it is None)
            # or overwrite_assignments is True, write the classification
            # assignment
            if current_classification is None or overwrite_assignments:
                classifications[key] = (major_classification,
                                        minor_classification)
            # Otherwise, write the major_classification only if the major
            # classification for the current_classification is the same
            elif current_classification[0] == major_classification:
                classifications[key] = (major_classification,
                                        minor_classification)
            # Otherwise do nothing


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


def _is_helix(residue):
    """Return true if the residue's Ramachandran angles fall within helical 
    values."""
    phi, psi = residue.ramachandran_angles
    if phi is None or psi is None:
        return False

    # Determine if the phi/psi angles are within a sheet range.
    phi_min, phi_max = settings.helix_phi
    psi_min, psi_max = settings.helix_psi

    return phi_min <= phi <= phi_max and psi_min <= psi <= psi_max