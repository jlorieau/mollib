"""
Functions to add hydrogens to atoms in a molecule.

Implementation
--------------
The current implementation has the following features:

- adds hydrogens in an idealized geometry with respect to local geometry
  groups.
- Ionizeable groups, like Asp and Glu carboxylates, are protonated based on the
  molecule's sample pH and the set pKa values. When more than one heavy atom
  can be protonated, the heavy atom closest to a hydrogen bond accector will
  be selected.
- Molecules with heteroatoms are correctly hydrogenated, provided the
  topological information is provided.

Written Molecule Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~
Category: 'Structural Features'
    hydrogenate: bool
        True if the molecule was hydrogenated

"""
# TODO: work with the residue ionization groups
# TODO: Add framework to modify ideal hydrogen geometries

import logging
from math import sqrt, cos, sin, pi
from itertools import chain

from mollib.core import calc_vector, vector_length, within_distance, cross
from . import settings


def add_hydrogens_to_residue(residue):
    """Add hydrogens to the given residue.

    Any ion groups for this residue are returned for processing last in the
    molecule.
    """
    # Get a list of the ionizeable groups to figure out how many protons
    # they need separately
    ion_groups = residue.ionizeable_groups
    ion_atoms = list(chain(*[i.possible_atoms for i in ion_groups]))

    for atom in residue.atoms:
        if atom.element == 'H' or atom.element == 'D':
            continue

        # Process ionizeable groups separately
        if atom in ion_atoms:
            continue

        return_value = add_hydrogen_to_atom(atom)
        if return_value is False:
            msg = "Could not add hydrogen to '{}' ".format(atom.fullname)
            logging.warning(msg)

    return ion_groups


def add_hydrogens_to_queue(q, result):
    """Add hydrogens to a queue of residues.

    This function is used for the concurrent implementation of
    add_hydrogens_to_residue. (see above)
    """
    while True:
        residue = q.get()

        if residue is None:
            q.task_done()
            break

        pH = residue.molecule.pH

        # Get a list of the ionizeable groups to figure out how many protons
        # they need separately
        ion_groups = residue.ionizeable_groups
        ion_atoms = list(chain(*[i.possible_atoms for i in ion_groups]))

        for atom in residue.atoms:
            if atom.element == 'H' or atom.element == 'D':
                continue

            # Process ionizeable groups separately
            if atom in ion_atoms:
                continue

            return_value = add_hydrogen_to_atom(atom)
            if return_value is False:
                msg = "Could not add hydrogen to '{}' ".format(atom.fullname)
                logging.warning(msg)

        # Now process ionizeable groups, based on the molecule's pH
        add_hydrogen_iongroups(pH=pH, ion_groups=ion_groups)
        q.task_done()
        result.put(residue)


def add_hydrogens(molecule, strip=True, force=False):
    """Add hydrogens to a molecule, if it hasn't been done already.

    Parameters
    ----------
    molecule: :obj:`mollib.Molecule`
        Molecule object to add hydrogens to.
    strip: bool, optional
        If True, all hydrogen atoms will be stripped from the molecule before
        adding new hydrogens
    force: bool, optional
        If True, the molecule will be rehydrogenated, even if it has before.
    """
    # TODO: Make parralel with dict residue, namedtuple atom
    # NOTE: I tried multiprocessing this function. It reaches the recursion
    # limit on large molecules, and the ionization doesn't work.

    already_hydrogenated = molecule.get_parameter('Structural Features',
                                                  'hydrogenate')
    if already_hydrogenated and not force:
        return None

    if strip:
        molecule.strip_atoms(element='H')

    pH = molecule.pH

    # Implementation #1 - single process
    ion_groups_list = []
    for residue in molecule.residues:
        ion_groups_list.append(add_hydrogens_to_residue(residue))

    # Process ion_groups separately
    pH = residue.molecule.pH
    for ion_groups in ion_groups_list:
        add_hydrogen_iongroups(pH=pH, ion_groups=ion_groups)

    # Now that all of the atoms have been added, clear the cache
    molecule.clear_cache(scope='add_atoms')

    # Mark that the molecule was hydrogenated
    molecule.set_parameter('Structural Features', 'hydrogenate', True)

    # Implementation #2 - multiprocess (twice as slow)
    #     The concurrency with map seems to work very well, but it chokes
    #     on residues within ionizeable groups, likely because the molecule
    #     has to be traversed for these.
    #
    # queue = mp.JoinableQueue()
    # result = mp.Queue()
    # num_workers = mp.cpu_count()
    #
    # workers = [mp.Process(target=add_hydrogens_to_queue, args=(queue, result))
    #            for i in range(num_workers)]
    # for w in workers:
    #     w.start()
    #
    # count = 0
    # for residue in molecule.residues:
    #     queue.put(residue)
    #     count +=1
    #
    # for i in range(num_workers):
    #     queue.put(None)
    #
    # queue.join()
    #
    # while count > 0:
    #     residue = result.get()
    #     molecule_residue = molecule[residue.chain.id][residue.number]
    #     for atom in residue.atoms:
    #         if atom.name not in molecule_residue:
    #             name = atom.name
    #             pos = atom.pos
    #             element = atom.element
    #             molecule.add_atom(name=name, pos=pos, element=element,
    #                               residue=molecule_residue)
    #     count -= 1

def add_hydrogen_to_atom(atom, number_hydrogens=None):
    """Add hydrogens to atom, making a decision on its hybridization.

    Parameters
    ----------
    atom: :obj:`atom`
        The atom object to add one or more hydrogens to.
    number_hydrogens: int (optional)
        If specified, the given number of hydrogens will be added.

    Returns
    -------
    bool
        True, if the addition was successful
        False, if the addition was not successful
    """
    topology = atom.topology
    number_hydrogens = number_hydrogens or len([i for i in topology
                                               if i.startswith('H')])
    number_heavy_atoms = len([i for i in topology
                              if not i.startswith('H')])

    if atom.element == 'O' and number_hydrogens == 1:
        add_one_sp3_h(atom, settings.bond_length['O-H'])
    elif atom.element == 'N':
        if number_hydrogens == 1:
            return add_one_sp2_h(atom, settings.bond_length['N-H'])
        if number_hydrogens == 2:
            return add_two_sp2_h(atom, settings.bond_length['N-H2'])
        if number_hydrogens == 3:
            return add_three_sp3_h(atom, settings.bond_length['N-H'])
    elif atom.element == 'C':
        if number_heavy_atoms <= 2 and number_hydrogens == 1:
            return add_one_sp2_h(atom, settings.bond_length['C-H'])
        if number_heavy_atoms == 3 and number_hydrogens == 1:
            return add_one_sp3_h(atom, settings.bond_length['C-H'])
        if number_heavy_atoms == 2 and number_hydrogens == 2:
            return add_two_sp3_h(atom, settings.bond_length['C-H'])
        if number_hydrogens == 3:
            return add_three_sp3_h(atom, settings.bond_length['C-H'])


def add_one_sp2_h(atom, bond_length):
    """Add a single hydrogens to an sp2 hybridized atom.

    This function is useful for adding protons to add protons with a
    120-degree geometry, like backbone HNs.

    Parameters
    ----------
    atom : :obj:`atom`
        The atom to add hydrogens to.
    bond_length: float
        The length of the atom-h bond (in Angstroms).

    Returns
    -------
    bool:
        True if atom was successfully added, False if it wasn't.

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
    bonded_heavy_atoms = atom.bonded_heavy_atoms(sorted=True)
    residue = atom.residue
    molecule = atom.molecule

    # Get the name of the hydrogen to add
    h_name = 'H' + atom.name[1:]

    # The current implementation determines the geometry based on two bonded
    # atoms.
    if len(bonded_heavy_atoms) == 2:
        # Calculate the v1, v2 and bisector vectors
        v1 = calc_vector(atom.pos, bonded_heavy_atoms[0].pos, normalize=True)
        v2 = calc_vector(atom.pos, bonded_heavy_atoms[1].pos, normalize=True)
        bisect = v1 + v2
        length = vector_length(bisect)
        bisect /= length

        normal_outplane = cross(v1,v2)
        normal_outplane /= vector_length(normal_outplane)

        normal_inplane = cross(normal_outplane, bisect)
        normal_inplane /= vector_length(normal_inplane)

        # calculate the h position along the bisector
        # cos_inplane = np.cos(in_plane * np.pi / 180.)
        # sin_inplane = np.sin(in_plane * np.pi / 180.)
        # cos_outplane = np.cos(out_plane * np.pi / 180.)
        # sin_outplane = np.sin(out_plane * np.pi / 180.)
        #
        # h = (cos_inplane * sin_outplane * normal_inplane +
        #      sin_inplane * sin_outplane * normal_outplane +
        #      cos_outplane * bisect)
        h = bisect * bond_length
        h += atom.pos

        # Create the new hydrogen atom
        molecule.add_atom(name=h_name, pos=h, element='H',
                          residue=residue, clear_cache=False)
        return True

    # If only one bonded_heavy_atom is available (like in a CO oxygen), the
    # position is inferred by looking at the the heavy atoms bonded to the
    # bonded atoms.
    elif len(bonded_heavy_atoms) == 1:
        bonded = bonded_heavy_atoms[0]

        # Find the atoms bonded to the bonded heavy atom that aren't 'atom'
        bonded_to_bonded = bonded.bonded_heavy_atoms(sorted=True)
        bonded_to_bonded = [a for a in bonded_to_bonded if a != atom]
        if len(bonded_to_bonded) < 1:
            return False
        bonded_to_bonded = bonded_to_bonded[0]

        # This algorithm will place the proton trans from the bonded_to_bonded
        # atom.
        sin60 = sqrt(3.)/2.
        cos60 = 0.5

        v1 = calc_vector(atom.pos, bonded_heavy_atoms[0].pos, normalize=True)
        v2 = calc_vector(bonded_heavy_atoms[0].pos, bonded_to_bonded.pos,
                         normalize=True)

        n_v1v2 = cross(v1, v2)
        n_v1v2 /= vector_length(n_v1v2)

        n = cross(v1, n_v1v2)
        n /= vector_length(n)

        h = (v1 * cos60 - n * sin60) * bond_length
        h += atom.pos

        # Create the new hydrogen atom
        molecule.add_atom(name=h_name, pos=h, element='H',
                          residue=residue, clear_cache=False)
        return True
    else:
        logging.warning("Number of bonded atoms "
                        "for {} is {}.".format(atom, len(bonded_heavy_atoms)))
        return False


def add_two_sp2_h(atom, bond_length, jbnmr_convention=True):
    """Add two hydrogens to an sp2 hybridized atom.

    Parameters
    ----------
    atom : :obj:`atom`
        The atom object to add a hydrogen to.
    bond_length: float
        The length of the bond (in Angstroms) between the new proton an
        target_atom
    jbnmr_convention: bool, optional
        Use the JBNMR 12,1 (1998) naming convention.

    Returns
    -------
    bool
        True if atom was successfully added, False if it wasn't.


    .. note:: Protons added to double-bonded atoms will respect the E/Z
              convention for E:H1 and Z: H2. The exception is for Arginines,
              for which this convention is reversed.

    Examples
    --------
    >>> from mollib.core import Molecule, measure_angle, measure_distance
    >>> mol = Molecule('2MJB')
    >>> mol.strip_atoms(element='H')
    >>> Q2 = mol['A'][2]
    >>> ne2 = Q2['NE2']
    >>> add_two_sp2_h(ne2, 1.0)
    True
    >>> h1, h2, cd, oe1 = Q2['HE21'], Q2['HE22'], Q2['CD'], Q2['OE1']
    >>> angle1 = measure_angle(h1, ne2, cd)
    >>> angle2 = measure_angle(h2, ne2, cd)
    >>> print("{:.1f} {:.1f} degs".format(angle1, angle2))
    120.0 120.0 degs
    >>> d1 = measure_distance(h1, oe1)
    >>> d2 = measure_distance(h2, oe1)
    >>> print("E:{:.1f} A, Z:{:.1f} A".format(d1, d2))  # stereoassn't
    E:3.1 A, Z:2.5 A
    """
    bonded_heavy_atoms = atom.bonded_heavy_atoms(sorted=True)
    residue = atom.residue
    molecule = atom.molecule

    # Get the name of the hydrogen to add
    h_name_1 = 'H' + atom.name[1:] + '1'  # E
    h_name_2 = 'H' + atom.name[1:] + '2'  # Z

    # We use H1 for E and H2 for Z. However, if the JBNMR convention
    # is followed, some of the atoms have to be switched.
    # The following atoms have reversed stereo-specific assignment
    res_name = atom.residue.name
    atom_name = atom.name
    reversed = False

    if jbnmr_convention:
        if ((res_name == 'ARG' and (atom_name == 'NH1' or atom_name == 'NH2'))
            ):
            reversed = True

    if len(bonded_heavy_atoms) == 1:
        bonded = bonded_heavy_atoms[0]

        # Find the atoms bonded to the bonded heavy atom that aren't 'atom'
        bonded_to_bonded = bonded.bonded_heavy_atoms(sorted=True)
        bonded_to_bonded = [a for a in bonded_to_bonded if a != atom]

        if len(bonded_to_bonded) < 1:
            return False
        bonded_to_bonded = bonded_to_bonded[0]

        # Calculate the atom--bonded_atom vector (x-axis)
        x = calc_vector(atom.pos, bonded.pos, normalize=True)

        # Calculate the bonded_atom--bonded_to_bonded vector. The
        # bonded_to_bonded should be sorted by mass in order to preserve
        # the E/Z notation.
        vec = calc_vector(bonded.pos, bonded_to_bonded.pos, normalize=True)

        # Calculate the z-axis as the orthogonal vector to these two
        z = cross(x, vec)
        z /= vector_length(z)

        # Calculate the y-axis as the orthogonal to the x- and y- axes.
        # The y-axis is pointing toward the bonded_to_bonded atom for
        # the 'Z' assignment
        y = cross(z, x)
        y /= vector_length(y)

        # The new protons are along the x- and y-axes, tilted by 60-degrees
        cos_60 = 0.5
        sin_60 = sqrt(3.) / 2.

        # Get the vectors for each hydrogen
        if not reversed:
            h1_vec = x * cos_60 + y * sin_60
            h2_vec = x * cos_60 - y * sin_60
        else:
            h1_vec = x * cos_60 - y * sin_60
            h2_vec = x * cos_60 + y * sin_60

        h1 = h1_vec * bond_length + atom.pos
        h2 = h2_vec * bond_length + atom.pos

        # Add the new protons to the target_atom
        molecule.add_atom(name=h_name_1, pos=h1,
                          element='H', residue=atom.residue, clear_cache=False)
        molecule.add_atom(name=h_name_2, pos=h2,
                          element='H', residue=atom.residue, clear_cache=False)

        return True


def add_one_sp3_h(atom, bond_length):
    """Add a single hydrogen to an sp3 hybridized atom.

    This function is useful for adding methine protons, like backbone HAs.

    Parameters
    ----------
    atom : :obj:`atom`
        The atom to add a hydrogen to.
    bond_length: float
        The length of the atom-h bond (in Angstroms).

    Returns
    -------
    bool:
        True if atom was successfully added, False if it wasn't.


    Examples
    --------
    >>> from mollib.core import Molecule, measure_angle
    >>> mol = Molecule('2KXA')
    >>> mol.strip_atoms(element='H')
    >>> F3 = mol['A'][3]
    >>> n, ca, c, cb = F3['N'], F3['CA'], F3['C'], F3['CB']
    >>> add_one_sp3_h(ca, 1.0)
    True
    >>> ha = F3['HA']
    >>> angle1 = measure_angle(n, ca, ha)
    >>> angle2 = measure_angle(c, ca, ha)
    >>> angle3 = measure_angle(cb, ca, ha)
    >>> print("{:.1f} {:.1f} {:.1f} degs".format(angle1, angle2, angle3))
    108.4 108.7 109.5 degs
    """
    bonded_heavy_atoms = atom.bonded_heavy_atoms(sorted=True)
    molecule = atom.molecule

    # Get the name of the hydrogen to add
    h_name = 'H' + atom.name[1:]

    if len(bonded_heavy_atoms) == 3:
        # This is the typical situation for a CA-HA, for example.

        # Calculate the plane and the plane normal formed by
        # atom_1, atom_2 and atom_3
        v1 = calc_vector(bonded_heavy_atoms[0].pos, atom.pos, normalize=True)
        v2 = calc_vector(bonded_heavy_atoms[1].pos, atom.pos, normalize=True)
        v3 = calc_vector(bonded_heavy_atoms[2].pos, atom.pos, normalize=True)
        h_v = v1+v2+v3
        h_v /= vector_length(h_v)

        # calculate the h position along the norm
        h = -1.*h_v * bond_length + atom.pos

        # Create the new hydrogen atom
        molecule.add_atom(name=h_name, pos=h, element='H',
                          residue=atom.residue, clear_cache=False)
        return True
    elif len(bonded_heavy_atoms) == 1:
        # This might occur for an O-H hydroxyl group. To get the correct
        # geometry, the bonded-to-bonded atoms are needed
        bonded = bonded_heavy_atoms[0]

        # Find the atoms bonded to the bonded heavy atom that aren't 'atom'
        bonded_to_bonded = bonded.bonded_heavy_atoms(sorted=True)
        bonded_to_bonded = [a for a in bonded_to_bonded if a != atom]

        if len(bonded_to_bonded) < 1:
            return False
        bonded_to_bonded = bonded_to_bonded[0]

        # By default, make the new H trans to the bonded--bonded_to_bonded,
        # but in the same plane.
        v1 = calc_vector(bonded.pos, atom.pos, normalize=True)
        v2 = calc_vector(bonded.pos, bonded_to_bonded.pos, normalize=True)

        z = cross(v1, v2)
        z /= vector_length(z)

        y = cross(v1, z)
        y /= vector_length(y)

        h = (cos(19.5*pi/180.)*y -
             sin(19.5*pi/180.)*v1)  # 19.5 = 109.5 - 90deg
        h *= bond_length
        h += atom.pos

        # Create the new hydrogen atom
        molecule.add_atom(name=h_name, pos=h, element='H',
                          residue=atom.residue, clear_cache=False)
        return True
    else:
        msg = 'add_one_sp3_h() requires 2 or 3 heavy atoms for atom {}.'
        msg = msg.format(atom)
        logging.warning(msg)


def add_two_sp3_h(atom, bond_length, jbnmr_convention=True):
    """Add two hydrogens to an sp3 hybridized atom.

    This function is useful for adding methylenes, like backbone HA2/HA3 for
    glycines. The atoms are added stereospecifically with HA2 as pro-R and
    HA3 as pro-S.

    Parameters
    ----------
    atom: :obj:`atom`
        The :obj:`atom` object to add two hydrogens to.
    bond_length: float
        The length of the bond (in Angstrom) between the new proton and
        atom.
    jbnmr_convention: bool, optional
        Use the JBNMR 12,1 (1998) naming convention.

    Returns
    -------
    bool:
        True if atom was successfully added, False if it wasn't.

    Examples
    --------
    >>> from mollib.core import Molecule, measure_angle
    >>> mol = Molecule('2KXA')
    >>> mol.strip_atoms(element='H')
    >>> G4 = mol['A'][4]
    >>> n, ca, c = G4['N'], G4['CA'], G4['C']
    >>> add_two_sp3_h(ca, 1.0)
    True
    >>> ha2, ha3 = G4['HA2'], G4['HA3']
    >>> angle1 = measure_angle(ha2, ca, ha3)
    >>> angle2 = measure_angle(n, ca, ha3)
    >>> print('{:.1f} {:.1f} degs'.format(angle1, angle2))
    109.5 109.1 degs
    """
    bonded_heavy_atoms = atom.bonded_heavy_atoms(sorted=True)
    molecule = atom.molecule

    # Get the name of the hydrogen to add
    h_name_2 = 'H' + atom.name[1:] + '2'
    h_name_3 = 'H' + atom.name[1:] + '3'

    # We use H2 for pro-R and H3 for pro-S. However, if the JBNMR convention
    # is followed, some of the atoms have to be switched.
    # The following atoms have reversed stereo-specific assignment
    res_name = atom.residue.name
    atom_name = atom.name
    reversed = False

    if jbnmr_convention:
        # Note that this implementation depends on the sort_atom_list
        # implementation, which currently makes some mistakes in CG/CD atoms.
        if ((res_name == 'MET' and (atom_name == 'CB' or atom_name == 'CG')) or
            (res_name == 'GLN' and atom_name == 'CG') or
            (res_name == 'GLU' and atom_name == 'CG') or
            (res_name == 'LYS' and atom_name == 'CE') or
            (res_name == 'PRO' and atom_name == 'CD') or
            (res_name == 'SER' and atom_name == 'CB') or
            (res_name == 'ASP' and atom_name == 'CB') or
            (res_name == 'ASN' and atom_name == 'CB') or
            (res_name == 'ARG' and atom_name == 'CD') or
            (res_name == 'CYS' and atom_name == 'CB')
            ):
            reversed = True

    if len(bonded_heavy_atoms) == 2:
        # Calculate the bisecting vector between bonded[0]--atom--bonded[1]
        v1 = calc_vector(atom.pos, bonded_heavy_atoms[0].pos)
        v2 = calc_vector(atom.pos, bonded_heavy_atoms[1].pos)
        bisect = v1 + v2
        bisect /= vector_length(bisect)

        # Calculate the vector orthogonal to the bonded[0]--atom--bonded[1]
        # plane
        orthog = cross(v1, bisect)
        orthog /= vector_length(orthog)

        # Find the unit vector for hydrogens 2 and 3.
        s = 0.8166415551616789  # sin(54.75deg) -- 54.75 deg = 109.5/2 deg
        c = 0.5771451900372336  # cos(54.75deg)

        if reversed:
            h_v2 = bisect * c + orthog * s
            h_v2 /= vector_length(h_v2)

            h_v3 = bisect * c - orthog * s
            h_v3 /= vector_length(h_v3)
        else:
            h_v2 = bisect * c - orthog * s
            h_v2 /= vector_length(h_v2)

            h_v3 = bisect * c + orthog * s
            h_v3 /= vector_length(h_v3)

        # The new protons point along h_v1 and h_v2 with a length of bond_length
        # from target_atom
        h2 = h_v2 * bond_length + atom.pos  # pro-R
        h3 = h_v3 * bond_length + atom.pos  # pro-S

        # Create the new hydrogen atoms
        molecule.add_atom(name=h_name_2, pos=h2, element='H',
                          residue=atom.residue, clear_cache=False)
        molecule.add_atom(name=h_name_3, pos=h3, element='H',
                          residue=atom.residue, clear_cache=False)
        return True

    else:
        return False


def add_three_sp3_h(atom, bond_length, alpha=None):
    """Add three hydrogens to an sp3 hybridized atom.

    This function is useful for adding methyls, like HBs of alanines.

    Parameters
    ----------
    atom: :obj:`atom`
        The :obj:`atom` object to add hydrogens to.
    bond_length: float
        The length of the bond between atom and the new hydrogens.
    alpha: float, optional
        The rotation angle to use for the atom-H3 group in degrees.
        If specified, this angle will be used. Otherwise, if a value exists in
        the :obj:`Molecule`'s parameters, that value will be used. The value
        defaults to 0.0 degrees.


    .. note: This function will replace the alpha value with the one stored in
             the molecule property if it is set.

    Returns
    -------
    bool:
        True if atom was successfully added, False if it wasn't.

    Examples
    --------
    >>> from mollib.core import Molecule, measure_angle, measure_distance
    >>> mol = Molecule('2KXA')
    >>> mol.strip_atoms(element='H')
    >>> A5 = mol['A'][5]
    >>> n, ca, cb = A5['N'], A5['CA'], A5['CB']
    >>> add_three_sp3_h(cb, 1.0)
    True
    >>> hb1, hb2, hb3 = A5['HB1'], A5['HB2'], A5['HB3']
    >>> angle1 = measure_angle(hb1, cb, ca)
    >>> angle2 = measure_angle(hb2, cb, ca)
    >>> angle3 = measure_angle(hb3, cb, ca)
    >>> print('{:.1f} {:.1f} {:.1f} degs'.format(angle1, angle2, angle3))
    109.5 109.5 109.5 degs
    >>> d1 = measure_distance(hb1, cb)
    >>> d2 = measure_distance(hb2, cb)
    >>> d3 = measure_distance(hb3, cb)
    >>> print('{:.2f} {:.2f} {:.2f} A'.format(d1, d2, d3))
    1.00 1.00 1.00 A
    """
    # TODO: Add molecule parameters using a more human-readable name
    # TODO: Add beta tetrahedral distortion angle.
    # Set the alpha angle value
    molecule = atom.molecule
    alpha_parameter = molecule.get_parameter('Add_hydrogens',
                                             atom.fullname + '_alpha')
    alpha = alpha or alpha_parameter or 0.0

    # Get the heavy atoms bonded to this atom
    bonded_heavy_atoms = atom.bonded_heavy_atoms(sorted=True)

    # Get the name of the hydrogen to add
    h_name_1 = 'H' + atom.name[1:] + '1'
    h_name_2 = 'H' + atom.name[1:] + '2'
    h_name_3 = 'H' + atom.name[1:] + '3'

    if len(bonded_heavy_atoms) == 1:
        bonded = bonded_heavy_atoms[0]

        # Find the atoms bonded to the bonded heavy atom that aren't 'atom'
        bonded_to_bonded = bonded.bonded_heavy_atoms(sorted=True)
        bonded_to_bonded = [a for a in bonded_to_bonded if a != atom]

        if len(bonded_to_bonded) < 1:
            return False
        bonded_to_bonded = bonded_to_bonded[0]

        # Calculate the atom--bonded vector
        v1 = calc_vector(atom.pos, bonded.pos)
        v1 /= vector_length(v1)

        # Calculate the bonded--bonded_to_bonded vector
        v2 = calc_vector(bonded.pos, bonded_to_bonded.pos)
        v2 /= vector_length(v2)

        # Calculate the normal to v1 and v2
        norm1 = cross(v1, v2)
        norm1 /= vector_length(norm1)

        # Calculate the normal to norm1 and v1. v1, norm1 and norm2 form an
        # orthonormal set of vectors
        norm2 = cross(norm1, v1)
        norm2 /= vector_length(norm2)

        # All three new protons are tilted away from v1. However, the first
        # is opposite to v2 to make a trans group
        c_705 = 0.33380  # cos(70.5-deg)
        s_705 = 0.94264  # sin(70.5-deg)

        h_v1 = c_705 * v1 + s_705 * (norm2 * cos((alpha + 0.)*pi/180.) +
                                     norm1 * sin((alpha + 0.)*pi/180.))
        h1 = h_v1 * bond_length + atom.pos

        h_v2 = c_705 * v1 + s_705 * (norm2 * cos((alpha + 240.)*pi/180.) +
                                     norm1 * sin((alpha + 240.)*pi/180.))
        h2 = h_v2 * bond_length + atom.pos

        h_v3 = c_705 * v1 + s_705 * (norm2 * cos((alpha + 120.)*pi/180.) +
                                     norm1 * sin((alpha + 120.)*pi/180.))
        h3 = h_v3 * bond_length + atom.pos

        # Create the new hydrogen atoms
        molecule.add_atom(name=h_name_1, pos=h1, element='H',
                          residue=atom.residue, clear_cache=False)
        molecule.add_atom(name=h_name_2, pos=h2, element='H',
                          residue=atom.residue, clear_cache=False)
        molecule.add_atom(name=h_name_3, pos=h3, element='H',
                          residue=atom.residue, clear_cache=False)
        return True

    else:
        return False


def add_hydrogen_iongroups(pH, ion_groups):
    """Add hydrogen and make choices on which hydrogen to add for the given
    ionization groups based on pH and pKs.

    Parameters
    ----------
    pH: float
        The pH of the sample the molecule is in.
    ion_groups: list
        The list of ion_group named tuples to add hydrogens to.
        The named tuple has 'possible_atoms', which is a list of possible
        atoms to add hydrogens to, and 'pKs', which is a tuple of pKs


    .. note:: This function protonates atoms that are closer to hydrogen bond
              acceptors ('O' and 'N' atoms) before protonating other sites.
    """
    # A dict with atoms (keys) and the number of hydrogens to add
    atom_no_hydrogens = {}

    for ion_group in ion_groups:
        possible_atoms = ion_group.possible_atoms
        pKs = ion_group.pKs

        # Figure out how many hydrogens need to be added. This corresponds to
        # the number of hydrogens that are below the pK.
        number_of_hydrogens = len([pK for pK in pKs if pK > pH])

        # Now determine which of the possible_atoms need to be protonated.
        # This decision is made based on which group is closest to an H-bond
        # acceptor, like an 'O' and an 'N'. The closest_NO is a dict with lists
        # of atoms within 5.0 Angstroms from other residues. The atoms in the
        # list are sorted in reverse order by distance. (i.e. the last item is
        # the closest)
        closest_NO = {}
        for atom in possible_atoms:
            sorted_by_distance = sorted(within_distance(atom=atom,
                                                    cutoff=5.0,
                                                    exclude_intraresidue=True),
                                        key = lambda t: t[1],
                                        reverse=True)  # sort by dist. (reverse)
            closest_NO[atom] = sorted_by_distance

        # Go through and protonate the atoms
        atom_no_hydrogens = {atom:0 for atom in closest_NO.keys()}
        for count in range(number_of_hydrogens):
            # Find the possible atom with the closest hbond acceptor
            atom_with_closest = None
            closest_distance = None

            for atom, dist_list in sorted(closest_NO.items()):
                # Set the new closest_distance if None has been set yet
                # or if a new closest distance has been found
                if (closest_distance is None or
                    (len(dist_list) > 1 and
                     dist_list[-1] < closest_distance)):
                    atom_with_closest = atom

            # atom_with_closest should be set at this point, either with the
            # first atom in the closest_NO dict (sorted) or with the actual
            # atom_with_closest hbond acceptor
            no_hydrogens = (atom_no_hydrogens.setdefault(atom_with_closest, 0)
                            + 1)
            atom_no_hydrogens[atom_with_closest] = no_hydrogens

        # Redistribute multiple hydrogens if there are many atoms. This
        # keeps sidechains like His from getting two hydrogens on one
        # heavy atom.
        # FIXME: This rebalancing block isn't very robust and may not work
        # for other types of molecules
        if (sum(atom_no_hydrogens.values()) > 1 and
            sum(atom_no_hydrogens.values()) == len(atom_no_hydrogens)):
            atom_no_hydrogens = {atom:1 for atom in atom_no_hydrogens}

        for atom, no_hydrogens in atom_no_hydrogens.items():
            if no_hydrogens == 0:
                continue
            msg = "Added {} hydrogens to the ionizeable {}."
            logging.info(msg.format(no_hydrogens, atom))
            add_hydrogen_to_atom(atom=atom, number_hydrogens=no_hydrogens)

