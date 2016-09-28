"""
Classify hydrogen bonds and residues
"""
from . import settings as s

from itertools import groupby


def within_range(value, value_range, wrap=None):
    """Test whether the value is within the range tuple.

    Parameters
    ----------
    value: float or in
        The value to test
    value_range: tuple
        The tuple for the min/max range of values to test within
    wrap: float, optional
        If specified, testing will also be done for value +/- wrap. This is
        useful for angles that are periodic by 180 deg or 360 degrees.


    .. note:: A value of None is possible--Ramachandran angles of the first
              of last residue, for example. In this case, this function will
              return False.

    Examples
    --------
    >>> within_range(2.5, (2.4, 3.6))
    True
    >>> within_range(2.5, (2.6, 3.6))
    False
    """
    if value is None:
        return False
    if wrap is not None:
        return any((value_range[0] <= value <= value_range[1],
                    value_range[0] <= (value + wrap) <= value_range[1],
                    value_range[0] <= (value - wrap) <= value_range[1],))
    return (value_range[0] <= value <= value_range[1])

def classify_major(hbond):
    """Mark the major classification of hbonds.

    Sets the 'major_classification' attribute of
    :obj:`mollib.hbonds.HydrogenBond` objects.

    Parameters
    ----------
    hbond: :obj:`mollib.hbonds.HydrogenBond`
        The hydrogen bond to classify.

    Returns
    -------
    bool
        True if the classification was succesfully found, False otherwise.
    """
    d1, d2 = hbond.donor.atom1, hbond.donor.atom2
    a1, a2 = hbond.acceptor.atom1, hbond.acceptor.atom2

    classification = 'major_'
    donor_type = ''
    if d2.name == 'N' and d1.name in ('H', 'H1', 'H2', 'H3'):
        classification += 'bb_'
        donor_type = 'amide'
    elif d2.name == 'CA' and d1.name in ('HA', 'HA2', 'HA3'):
        classification += 'bb_'
        donor_type = 'aliphatic'
    elif d1.element == 'H' and d2.element == 'N':
        classification += 'sc_'
        donor_type = 'amide'
    elif d1.element == 'H' and d2.element == 'C':
        classification += 'sc_'
        donor_type = 'aliphatic'
    elif d1.element == 'H' and d2.element == 'O':
        classification += 'sc_'
        donor_type = 'hydroxyl'

    if a1.name == 'O' and a2.name == 'C':
        classification += 'bb_'
    elif a1.name == 'O' and a2.name == 'N':
        classification += 'bb_'
    else:
        classification += 'sc_'

    classification += donor_type

    if hasattr(s, classification):
        hbond.major_classification = getattr(s, classification)
        return True
    else:
        return False

def classify_minor(hbond):
    """Mark the minor classification of hbonds.

    Sets the 'minor_classification' attribute of
    :obj:`mollib.hbonds.HydrogenBond` objects. This function is run after
    classify_major.

    Parameters
    ----------
    hbond: :obj:`mollib.hbonds.HydrogenBond`
        The hydrogen bond to classify.

    Returns
    -------
    bool
        True if the classification was succesfully found, False otherwise.
    """
    try:
        donor_res = hbond.donor.atom2.residue
        acceptor_res = hbond.acceptor.atom2.residue
    except AttributeError:
        return False

    # All hydrogen bonds that are not Backbone-backbone amide are just
    # isolated hydrogen bonds.
    if hbond.major_classification != s.major_bb_bb_amide:
        hbond.minor_classification = s.minor_isolated
        return True

    # Calculate the difference in residue number for the H-bond donor
    # and acceptor
    delta = (donor_res.number - acceptor_res.number)

    # Get the residue i, i+1, i+2, i+3,
    res_i = acceptor_res
    res_i1 = (res_i.next_residue if res_i is not None else None)
    res_i2 = (res_i1.next_residue if res_i1 is not None else None)
    res_i3 = (res_i2.next_residue if res_i2 is not None else None)
    res_i4 = (res_i3.next_residue if res_i3 is not None else None)
    res_i5 = (res_i4.next_residue if res_i4 is not None else None)

    # Delta=3. These could be a beta-turn I, II or 310 helix, depending on
    # dihedrals.
    # Check that the residue i, i+1, i+2 and i+3 are all real residues
    # and not None
    if delta == 3 and all((res_i, res_i1, res_i2, res_i3)):
        # Get the Ramachandran phi/psi angles. The ramachandran_angles
        # method returns a tuple for (phi, psi)
        phi_psi = [r.ramachandran_angles for r in (res_i, res_i1, res_i2,
                                                   res_i3)]
        phi = [a[0] for a in phi_psi]
        psi = [a[1] for a in phi_psi]

        # Check the 310-helix
        if (all(within_range(a, s.helix_phi) for a in phi) and
            all(within_range(a, s.helix_psi) for a in psi)):
            hbond.minor_classification = s.minor_310
            return True

        # Check the Beta-turns. This part uses unions of sets to find the
        # common beta turn type for all of the dihedral angles.
        turn_type = ({k for k,v in s.beta_turn_i1_phi.items()
                      if within_range(phi[1], v)} &
                     {k for k,v in s.beta_turn_i1_psi.items()
                      if within_range(psi[1], v)} &
                     {k for k, v in s.beta_turn_i2_phi.items()
                      if within_range(phi[2], v)} &
                     {k for k, v in s.beta_turn_i2_psi.items()
                      if within_range(psi[2], v)})

        if len(turn_type) == 1:
            turn_type = turn_type.pop()
            classification = 'minor_beta_' + turn_type

            if hasattr(s, classification):
                hbond.minor_classification = getattr(s, classification)
                return True

    # Delta=4. These could be alpha-helix
    # Check that the residue i, i+1, i+2, i+3 and i+4 are all real residues
    # and not None
    elif delta==4 and all((res_i, res_i1, res_i2, res_i3, res_i4)):
        # Get the Ramachandran phi/psi angles. The ramachandran_angles
        # method returns a tuple for (phi, psi)
        phi_psi = [r.ramachandran_angles for r in (res_i, res_i1, res_i2,
                                                   res_i3, res_i4)]
        phi = [a[0] for a in phi_psi]
        psi = [a[1] for a in phi_psi]

        # Check the 310-helix
        if (all(within_range(a, s.helix_phi) for a in phi) and
                all(within_range(a, s.helix_psi) for a in psi)):
            hbond.minor_classification = s.minor_alpha
            return True

    # Delta=5. These could be pi-helix
    # Check that the residue i, i+1, i+2, i+3, i+4 and i+5 are all real
    # residues and not None
    elif delta==5 and all((res_i, res_i1, res_i2, res_i3, res_i4, res_i5)):

        # Get the Ramachandran phi/psi angles. The ramachandran_angles
        # method returns a tuple for (phi, psi)
        phi_psi = [r.ramachandran_angles for r in (res_i, res_i1, res_i2,
                                                   res_i3, res_i4, res_i5)]
        phi = [a[0] for a in phi_psi]
        psi = [a[1] for a in phi_psi]

        # Check the 310-helix
        if (all(within_range(a, s.helix_phi) for a in phi) and
                all(within_range(a, s.helix_psi) for a in psi)):
            hbond.minor_classification = s.minor_pi
            return True

    # At this point, helices and beta turns have been filtered out.
    # The only remaining options are parallel and anti-parallel beta-sheet
    # or isolated hydrogen bonds. The descrimination is based on torsion
    # angles
    d_phi, d_psi = donor_res.ramachandran_angles
    a_phi, a_psi = acceptor_res.ramachandran_angles

    # See of backbone dihedrals are consistent with a beta sheet
    if (abs(delta) >= 4 and   ## beta sheet minimum iter-residue spacing
        all(within_range(a, s.beta_phi, wrap=360.)
            for a in (d_phi, a_phi)) and
        all(within_range(a, s.beta_psi, wrap=360.)
            for a in (d_psi, a_psi))):

        # This is a beta sheet. The check for parallel or anti-parallel
        # is done elsewhere.
        hbond.minor_classification = s.minor_beta
        return True

    hbond.minor_classification = s.minor_isolated
    return True


def group_hbonds(hbonds):
    """Group hydrogen bonds into lists of list of hbonds with consecutive
    residue numbers.

    Parameters
    ----------
    hbonds:  hbonds: list of :obj:`mollib.hbonds.HydrogenBond`
        Hydrogen bonds to classify.

    Returns
    -------
    list of lists of :obj:`mollib.hbonds.HydrogenBond` objects.
    """
    key = lambda i: i[0] - i[1].acceptor.atom2.residue.number

    return [[i[1] for i in g] for k, g in
            groupby(enumerate(hbonds), key)]


def classify_groups(hbonds):
    """Classify groups of hbonds together.

    This function is run after classify_minor to further quality the
    minor classification in terms of parallel/antiparallel sheets, and the
    terminii of helices. Sets the 'minor_classification' attribute of
    :obj:`mollib.hbonds.HydrogenBond` objects.

    Parameters
    ----------
    hbonds: list of :obj:`mollib.hbonds.HydrogenBond`
        Hydrogen bonds to classify.

    Returns
    -------
    bool
        True if the classification was succesfully found, False otherwise.
    """
    # Sort the hbonds by residue number
    hbonds_sorted = sorted(hbonds,
                           key=lambda k: k.acceptor.atom2.residue.number)

    # filter the alpha-helices
    ahelix_hbonds = [hbond for hbond in hbonds
                     if hbond.minor_classification == s.minor_alpha]

    # Group the hbonds by contiguous stretches in residue number.
    ahelix_groups = group_hbonds(ahelix_hbonds)

    # Go through the contiguous stretches of hbonds and identify secondary
    # structure units.
    for hbonds_group in ahelix_groups:

        # Mark the starting and ending residues of the helices
        if len(hbonds_group) > 3:
            for hbond in hbonds_group[0:2]:
                hbond.minor_classification = s.minor_alpha_N
            for hbond in hbonds_group[-2:]:
                hbond.minor_classification = s.minor_alpha_C

    # Identify parallel and anti-parallel sheets
    # filter the beta-sheets
    # beta_hbonds = [hbond for hbond in hbonds
    #                  if hbond.minor_classification == s.minor_beta]
    #
    #

    return True


def classify_hbonds(hbonds):
    """Classify hydrogen bonds.
    """
    for hbond in hbonds:
        classify_major(hbond)
        classify_minor(hbond)
    classify_groups(hbonds)

