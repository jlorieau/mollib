"""
Classify hydrogen bonds
"""
from . import settings


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


class HbondClassifier(object):
    """Hydrogen bond classification class.

    The classifications are made with the classify_major_* and classify_minor_*
    methods, which return True if the hbond was classified, or False if it
    wasn't. The order of the functions are from most specific to most general.


    .. note:: This class is a singleton and will return the same instance on
              construction.
    """

    major_bb_bb_amide = 'bb-bb amide'
    major_bb_sc_amide = 'bb-sc amide'
    major_sc_bb_amide = 'sc-bb amide'
    major_sc_sc_amide = 'sc-sc amide'

    major_bb_bb_aliphatic = 'bb-bb aliph.'
    major_bb_sc_aliphatic = 'bb-sc aliph.'
    major_sc_bb_aliphatic = 'sc-bb aliph.'
    major_sc_sc_aliphatic = 'sc-sc aliph.'

    major_bb_sc_hydroxyl = 'bb-sc hydroxyl'
    major_sc_bb_hydroxyl = 'sc-bb hydroxyl'
    major_sc_sc_hydroxyl = 'sc-sc hydroxyl'

    minor_beta_turnI = "(type I turn)"
    minor_beta_turnII = "(type II turn)"
    minor_beta_turnIp = "(type I' turn)"
    minor_beta_turnIIp = "(type II' turn)"

    minor_beta_anti = '(sheet, anti-paralle;)'
    minor_beta_par = '(sheet, parallel)'

    minor_310 = '(310-helix)'
    minor_alpha = '(alpha-helix)'
    minor_pi = '(pi-helix)'
    minor_isolated = '(isolated)'

    _instance = None

    def __new__(cls, *args, **kwargs):
        "Create and return a singleton"
        if not cls._instance:
            cls._instance = super(HbondClassifier, cls).__new__(cls, *args,
                                                                **kwargs)
        return cls._instance

    def classify_major(self, hbond):
        """Mark the major classification of hbonds.

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
        if d1.name == 'H' and d2.name == 'N':
            classification += 'bb_'
            donor_type = 'amide'
        elif d1.name == 'HA' and d2.name == 'CA':
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

        if hasattr(self, classification):
            hbond.major_classification = getattr(self, classification)
            return True
        else:
            return False

    def classify_minor(self, hbond):
        """Mark the minor classification of hbonds.

        This function is run after classify_major.

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

            # Get the residu i, i+1, i+2, i+3, i+4 and i+5
            res_i = acceptor_res
            res_i1 = res_i.next_residue
            res_i2 = res_i1.next_residue
            res_i3 = res_i2.next_residue
            res_i4 = res_i3.next_residue
            res_i5 = res_i4.next_residue

        except AttributeError:
            return False

        # All hydrogen bonds that are not Backbone-backbone amide are just
        # isolated hydrogen bonds.
        if hbond.major_classification != self.major_bb_bb_amide:
            hbond.minor_classification = self.minor_isolated
            return True

        # Calculate the difference in residue number for the H-bond donor
        # and acceptor
        delta = (donor_res.number - acceptor_res.number)

        # Delta=3. These could be a beta-turn I, II or 310 helix, depending on
        # dihedrals.
        if delta == 3:
            # Get the Ramachandran phi/psi angles. The ramachandran_angles
            # method returns a tuple for (phi, psi)
            phi_psi = [r.ramachandran_angles for r in (res_i, res_i1, res_i2,
                                                       res_i3)]
            phi = [a[0] for a in phi_psi]
            psi = [a[1] for a in phi_psi]

            # Check the 310-helix
            if (all(within_range(a, settings.helix_phi) for a in phi) and
                all(within_range(a, settings.helix_psi) for a in psi)):
                hbond.minor_classification = self.minor_310
                return True

            # Check the Beta-turns. This part uses unions of sets to find the
            # common beta turn type for all of the dihedral angles.
            turn_type = ({k for k,v in settings.beta_turn_i1_phi.items()
                          if within_range(phi[1], v)} &
                         {k for k,v in settings.beta_turn_i1_psi.items()
                          if within_range(psi[1], v)} &
                         {k for k, v in settings.beta_turn_i2_phi.items()
                          if within_range(phi[2], v)} &
                         {k for k, v in settings.beta_turn_i2_psi.items()
                          if within_range(psi[2], v)})
            
            if len(turn_type) == 1:
                turn_type = turn_type.pop()
                classification = 'minor_beta_' + turn_type

                if hasattr(self, classification):
                    hbond.minor_classification = getattr(self, classification)
                    return True
                return False

        # Delta=4. These could be alpha-helix
        elif delta==4:
            # Get the Ramachandran phi/psi angles. The ramachandran_angles
            # method returns a tuple for (phi, psi)
            phi_psi = [r.ramachandran_angles for r in (res_i, res_i1, res_i2,
                                                       res_i3, res_i4)]
            phi = [a[0] for a in phi_psi]
            psi = [a[1] for a in phi_psi]

            # Check the 310-helix
            if (all(within_range(a, settings.helix_phi) for a in phi) and
                    all(within_range(a, settings.helix_psi) for a in psi)):
                hbond.minor_classification = self.minor_alpha
                return True

        # Delta=5. These could be pi-helix
        elif delta==5:
            # Get the Ramachandran phi/psi angles. The ramachandran_angles
            # method returns a tuple for (phi, psi)
            phi_psi = [r.ramachandran_angles for r in (res_i, res_i1, res_i2,
                                                       res_i3, res_i4, res_i5)]
            phi = [a[0] for a in phi_psi]
            psi = [a[1] for a in phi_psi]

            # Check the 310-helix
            if (all(within_range(a, settings.helix_phi) for a in phi) and
                    all(within_range(a, settings.helix_psi) for a in psi)):
                hbond.minor_classification = self.minor_pi
                return True

        # At this point, helices and beta turns have been filtered out.
        # The only remaining options are parallel and anti-parallel beta-sheet
        # or isolated hydrogen bonds. The descrimination is based on torsion
        # angles
        d_phi, d_psi = donor_res.ramachandran_angles
        a_phi, a_psi = acceptor_res.ramachandran_angles

        # See of backbone dihedrals are consistent with a beta sheet
        if (all(within_range(a, settings.beta_phi, wrap=360.)
                for a in (d_phi, a_phi)) and
            all(within_range(a, settings.beta_psi, wrap=360.)
                for a in (d_psi, a_psi))):

            # This is a beta sheet. Is it parallel or anti-parallel. We must
            # look up the residue numbers for the preceeding and proceeding
            # residues of both.
            try:
                acceptor_delta = (acceptor_res.next_residue.number -
                                  acceptor_res.prev_residue.number)
                donor_delta = (donor_res.next_residue.number -
                               donor_res.prev_residue.number)
            except AttributeError:
                return False

            # If acceptor_delta and donor_delta have the same sign, then
            # they're parallel. Otherwise they're anti-parallel
            if acceptor_delta * donor_delta > 0.:
                hbond.minor_classification = self.minor_beta_par
            else:
                hbond.minor_classification = self.minor_beta_anti
            return True

        hbond.minor_classification = self.minor_isolated
        return True

classifier = HbondClassifier()

def classify(hbond):
    """Classify the hydrogen bond.

    Returns
    -------
    bool
        True if the classification was succesful, otherwise False.
    """
    global classifier

    return_value = False
    return_value &= classifier.classify_major(hbond)
    return_value &= classifier.classify_minor(hbond)

    return return_value