"""
Functions to add hydrogens to atoms in a molecule.
"""
import logging

from mollib.core import calc_vector, vector_length


def add_hydrogens(molecule, strip=True):
    """

    Parameters
    ----------
    molecule
    strip

    Returns
    -------

    """
    raise NotImplementedError


def add_one_sp2_h(atom, bond_length):
    """Add a single hydrogens to an sp2 hybridized atom.

    This function is useful for adding protons to add protons with a
    120-degree geometry, like backbone HNs.

    Parameters
    ----------
    atom : :obj:`atom`
        The atom to add hydrogens to.

    Returns
    -------
    bool:
        True if atom was successfully added, False if it wasn't.


    .. note:: Protons added to double-bonded atoms will respect the (E), (Z)
              assignment convention.

    Examples
    --------
    >>> from mollib.core import Molecule, measure_angle
    >>> mol = Molecule('2KXA')
    >>> mol.strip_atoms(element='H')
    >>> F3 = mol['A'][3]
    >>> add_one_sp2_h(F3['N'], 1.0)
    True
    >>> hn = F3['H']
    >>> c_prev, n, ca = mol['A'][2]['C'], F3['N'], F3['CA']
    >>> angle1 = measure_angle(c_prev, n, hn)
    >>> angle2 = measure_angle(ca, n, hn)
    >>> print("{:.1f} {:.1f} degs".format(angle1, angle2))
    119.4 119.4 degs
    """
    bonded_atoms = [a for a in atom.bonded_atoms
                    if a.element != 'H' or a.element != 'D']
    bonded_atoms = sorted(bonded_atoms, key=lambda a: a.mass,
                          reverse=True)
    residue = atom.residue
    molecule = atom.molecule

    # Get the name of the hydrogen to add
    h_name = list(filter(lambda x:x.startswith('H'), atom.topology))
    if not h_name:  # Not eligible Hs found
        msg = ("Could not find a H atom for {}. ".format(atom) +
               "Eligible atoms are: {}".format(atom.bonded_atoms))
        logging.warning(msg)
        return False
    else:
        h_name = h_name[0]  # return the first hydrogen found

    # The current implementation determines the geometry based on two bonded
    # atoms. If only one is available (like in a CO oxygen), the position
    # is inferred by looking at the the heavy atoms bonded to the bonded atoms.
    if len(bonded_atoms) == 2:
        # Calculate the v1, v2 and bisector vectors
        v1 = calc_vector(atom, bonded_atoms[0])
        v2 = calc_vector(atom, bonded_atoms[1])
        bisect = v1 + v2
        length = vector_length(bisect)
        bisect /= length

        # calculate the h position along the bisector
        h = bisect * bond_length + atom.pos

        # Create the new hydrogen atom

        molecule.add_atom(name=h_name, pos=h, charge=0.0, element='H',
                          residue=residue)
        return True
    else:
        logging.warning("Number of bonded atoms "
                        "for {} is {}.".format(atom, len(bonded_atoms)))
        return False





