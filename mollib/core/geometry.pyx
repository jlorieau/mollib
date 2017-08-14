"""
Tools to measure geometries in molecules.
"""
# TODO: def measure_rmsd(molecule1, molecule2, atoms=None)

# cython: profile=False

from libc.math cimport sqrt, ceil, acos, atan2, M_PI as pi

import numpy as np
cimport numpy as np
cimport cython

from .geometry_box cimport Box
from .utils import filter_atoms


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef double vector_length(np.ndarray[np.float64_t, ndim=1] vector):
    """Return the length (in A) of a vector"""
    cdef Py_ssize_t i
    cdef Py_ssize_t vector_size = vector.shape[0]
    cdef double v2 = 0.0
    for i in range(vector_size):
        v2 += vector[i]*vector[i]
    return sqrt(v2)



@cython.boundscheck(False)
@cython.wraparound(False)
cpdef np.ndarray[np.float64_t, ndim=1] cross(np.ndarray[np.float64_t, ndim=1] a,
                                             np.ndarray[np.float64_t, ndim=1] b):
    """Return the cross product between two vectors."""
    cdef np.ndarray[np.float64_t, ndim=1] c = np.zeros((3), dtype=np.float64)
    c[0] = a[1]*b[2] - a[2]*b[1]
    c[1] = a[2]*b[0] - a[0]*b[2]
    c[2] = a[0]*b[1] - a[1]*b[0]
    return c


def calc_vector(vector_i, vector_j, normalize=True):
    """Returns the vector between atoms 'i' and 'j' with optional
    normalization."""
    vec = vector_i - vector_j

    if normalize:
        length = vector_length(vec)
        vec /= length
    return vec


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef double measure_distance(object atom_1, object atom_2):
    """Measure the atom_1--atom_2 distance.

    Parameters
    ----------
    atom_1 : :obj:`Atom`
        The first atom.
    atom_2 : :obj:`Atom`
        The second atom.

    Returns
    -------
    distance : float
        The distance (in Angstroms)

    Examples
    --------
    >>> from mollib.core import Molecule, measure_distance
    >>> mol = Molecule('2KXA')
    >>> d = measure_distance(mol['A'][3]['CA'], mol['A'][3]['HA'])
    >>> print("{:.2f} A".format(d))
    1.08 A

    """
    cdef double x, y, z
    cdef double length = 0.0
    cdef np.ndarray[np.float64_t, ndim=1] vi, vj

    vi = atom_1.pos
    vj = atom_2.pos

    x = vi[0] - vj[0]
    y = vi[1] - vj[1]
    z = vi[2] - vj[2]

    length = sqrt(x * x + y * y + z * z)
    return length

@cython.boundscheck(False)
@cython.wraparound(False)
cdef double _within_distance(double[:] point1, double[:] point2,
                           double distance_cutoff) nogil:
    """Test whether two points (specified by vectors) are within
    distance_cutoff.

    Parameters
    ----------
    point1: double[:]
        A 1D memoryview of the first point.
    point2
        A 2D memoryview of the second point
    Returns
    -------
    double
        - The distance between the two points if within the distance_cutoff
        - -1.0 if the distance between the two points is outside the
          distance_cutoff
    """
    # assert
    cdef double vx, vy, vz, d2
    vx = point1[0] - point2[0]
    vy = point1[1] - point2[1]
    vz = point1[2] - point2[2]

    d2 = vx * vx
    d2 += vy * vy
    d2 += vz * vz

    if d2 < distance_cutoff * distance_cutoff:
        return sqrt(d2)
    else:
        return -1.0


def within_distance(atom, cutoff, elements='', exclude_intraresidue=False):
    """Find all atoms of element within the specified distance (in Angstroms)
       of atom.

   Parameters
   ----------
   atom: :obj:`atom`
       The atom to find atoms around it.
   cutoff: float
       The distance boundary between atom and atoms of element to return.
   elements: str
       The element names of the atoms to return. This string supports the
       or character '|'.
       If '', all atoms within the distance will be returned
       ex: 'H|C|N' for all H, C and N atoms
   exclude_intraresidue: bool
       If True, atoms within the same residue as atom will be excluded.
   atom_selelction: iterable, optional
       If specified, the nearest neighbors will be searched from this iterable
       instead of the atom.molecule attribute.


       ..note : The elements and exclude_intraresidue methods are ignored when
                an atom_selection is specified.

   Returns
   -------
   list of tuples
       A list of tuples with (atom objects, distance).

    Examples
    --------
    >>> from mollib.core import Molecule, within_distance
    >>> mol = Molecule('2KXA')
    >>> within_distance(mol['A'][5]['CA'], cutoff=2.0)
    [A.A5.CB, A.A5.C, A.A5.N, A.A5.HA]

    """

    if 'box' not in atom.molecule.cache:
        box = Box(atom.molecule.atoms)
        atom.molecule.cache['box'] = box

    box = atom.molecule.cache['box']
    atoms = box.get_points(atom, cutoff)

    # Filter the atoms in the iterable
    element_list = elements.split('|') if elements != '' else []

    atoms = [a for a in atoms if
             ((a != atom) and
              (not element_list or a.element in element_list) and
              (not exclude_intraresidue or a.residue != atom.residue))]
    return atoms


def measure_angle(atom_1, atom_2, atom_3):
    """Measure the atom_1--atom_2--atom_3 angle.

    Parameters
    ----------
    atom_1: :obj:`Atom`
        The first atom of the angle to measure.
    atom_2: :obj:`Atom`
        The second atom of the angle to measure.
    atom_3: :obj:`Atom`
        The third atom of the angle to measure.

    Returns
    -------
    angle : float
        The angle (in degrees).

    Examples
    --------
    >>> from mollib.core import Molecule, measure_angle
    >>> mol = Molecule('2KXA')
    >>> gly4 = mol['A'][4]
    >>> angle = measure_angle(gly4['HA2'], gly4['CA'], gly4['HA3'])
    >>> print("{:.1f} deg".format(angle))
    109.3 deg

    """
    v1 = calc_vector(atom_2.pos, atom_1.pos, normalize=True)
    v2 = calc_vector(atom_2.pos, atom_3.pos, normalize=True)

    angle = acos(np.dot(v1, v2))
    return angle*180./pi


def measure_dihedral(atom_1, atom_2, atom_3, atom_4):
    """Measure the atom_1--atom_2--atom_3--atom_4 dihedral angle.

    Parameters
    ----------
    atom_1 : :obj:`Atom`
        The first atom.
    atom_2 : :obj:`Atom`
        The second atom.
    atom_3 : :obj:`Atom`
        The third atom.
    atom_4 : :obj:`Atom`
        The fourth atom.

    Returns
    -------
    angle : float
        The dihedral angle.

    Examples
    --------
    >>> from mollib.core import Molecule, measure_dihedral
    >>> mol = Molecule('2KXA')
    >>> F3 = mol['A'][3]
    >>> angle = measure_dihedral(F3['N'], F3['CA'], F3['CB'], F3['CG'])
    >>> print("{:.1f} deg".format(angle))
    -62.0 deg

    """
    # Calculate the normalized vectors.
    ab = calc_vector(atom_1.pos, atom_2.pos, normalize=True)
    bc = calc_vector(atom_2.pos, atom_3.pos, normalize=True)
    cd = calc_vector(atom_3.pos, atom_4.pos, normalize=True)

    # The dihedral is the angle between the plans a-b-c and b-c-d
    # The angle between these plans can be calculated from their
    # normals (cross products)
    cdef np.ndarray[np.float64_t, ndim=1] n1 = cross(ab, bc)
    cdef np.ndarray[np.float64_t, ndim=1] n2 = cross(bc, cd)

    # The angle between n1 and n2 can be calculated with acos. However
    # the followin atan2 relationship returns a number between 0 and
    # 2pi
    cdef np.ndarray[np.float64_t, ndim=1] m1 = cross(n1, bc)
    x = np.dot(n1, n2)
    y = np.dot(m1, n2)
    angle = atan2(y, x) * 180. / np.pi
    return angle


def measure_distances(molecule, selector1, selector2,
                      **filters):
    """Measure the distances for atoms selected by selector1 and selector2.

    Parameters
    ----------
    molecule: :obj:`mollib.Molecule`
        The molecule object from which distances will be measured.
    selector1: str
        A string for locator of the first atom(s) to measure the distances
        from.
    selector2: str
        A string for the locator of the second atom(s) to measure the distances
        to.
    filters: dict
        See :func:`mollib.core.filter_atoms`


        .. note:: The selectors respects the rules outlined in
                  :class:`mollib.Molecule.get_atoms`

    Returns
    -------
    list of tuples
        A list of tuples for (atom1, atom2, distance) for each distance
        selected.

    Examples
    --------
    >>> from mollib.core import Molecule, measure_distances
    >>> mol = Molecule('2MUV')
    >>> kwargs = {'residue_delta': 0, 'only_intra': True}
    >>> dists = measure_distances(mol, '23:25.N', '23:25.CA', **kwargs)
    >>> for i in dists: print(i)
    (A.S23.N, A.S23.CA, 1.46)
    (A.D24.N, A.D24.CA, 1.45)
    (A.P25.N, A.P25.CA, 1.49)
    >>> dists = measure_distances(mol, 'A:C.23.N', 'A:C.23.N')
    >>> for i in dists: print(i)
    (A.S23.N, B.S23.N, 11.23)
    (A.S23.N, C.S23.N, 15.49)
    (B.S23.N, C.S23.N, 11.04)
    >>> dists = measure_distances(mol, '23.N', '24.N', exclude_intra=True)
    >>> for i in dists: print(i)
    (A.S23.N, A.D24.N, 3.45)

    """
    # Returned dict
    results = {}

    # This function logs an error if a1 or a2 isn't properly
    # formatted. An additional message is not needed. Just skip
    # it if both atoms aren't found.
    atoms1 = molecule.get_atoms(selector1)
    atoms2 = molecule.get_atoms(selector2)

    if not atoms1 or not atoms2:
        return []

    for i in atoms1:
        for j in atoms2:
            # Sort the atoms into a tuple
            key = (i, j)

            # The atoms in forward or reverse order are duplicates
            if key in results or key[::-1] in results:
                continue

            # Process the filters
            if filter_atoms(*key, **filters):
                continue

            # Add the distance to the results
            dist = round(measure_distance(i, j), 2)
            results[key] = dist

    # Return the sorted list of tuples. Sort first by chain number, then
    # by residue number
    sort_key = lambda i: (i[0][0].chain.id,
                          i[0][1].chain.id,
                          i[0][0].residue.number,
                          i[0][1].residue.number)
    return [(k[0], k[1], v)
            for k, v in sorted(results.items(), key=sort_key)]


def measure_angles(molecule, selector1, selector2, selector3,
                   **filters):
    """Measure the angles for atoms selected by selector1, selector 2 and
    selector3.

    Parameters
    ----------
    molecule: :obj:`mollib.Molecule`
        The molecule object from which distances will be measured.
    selector1: str
        A string for the selector of the first atom(s) of the angle.
    selector2: str
        A string for the selector of the second atom(s) of the angle.
    selector3: str
        A string for the selector of the third atom(s) of the angle.
    filters: dict
        See :func:`mollib.core.filter_atoms`


        .. note:: The locator respects the rules outlined in
                  :class:`mollib.Molecule.get_atoms`

    Returns
    -------
    list of tuples
        A sorted list of tuples with (atom1, atom2, atom3, angle) for each
        angle selected.

    Examples
    --------
    >>> from mollib.core import Molecule, measure_angles
    >>> mol = Molecule('2KXA')
    >>> ang = measure_angles(mol, '3:5.N', '3:5.CA', '3:5.C', only_intra=True)
    >>> for i in ang: print(i)
    (A.F3.N, A.F3.CA, A.F3.C, 110.7)
    (A.G4.N, A.G4.CA, A.G4.C, 110.9)
    (A.A5.N, A.A5.CA, A.A5.C, 110.7)
    >>> kwargs = {'residue_delta': 1, 'bonded': True}
    >>> ang = measure_angles(mol, '3:8.C', '3:8.N', '3:8.CA', **kwargs)
    >>> for i in ang: print(i)
    (A.F3.C, A.G4.N, A.G4.CA, 120.8)
    (A.G4.C, A.A5.N, A.A5.CA, 121.5)
    (A.A5.C, A.I6.N, A.I6.CA, 120.7)
    (A.I6.C, A.A7.N, A.A7.CA, 121.2)
    (A.A7.C, A.G8.N, A.G8.CA, 120.9)

    """
    # Returned dict
    results = {}

    # This function logs an error if a1 or a2 isn't properly
    # formatted. An additional message is not needed. Just skip
    # it if both atoms aren't found.
    atoms1 = molecule.get_atoms(selector1)
    atoms2 = molecule.get_atoms(selector2)
    atoms3 = molecule.get_atoms(selector3)

    if not atoms1 or not atoms2 or not atoms3:
        return []

    for i in atoms1:
        for j in atoms2:
            for k in atoms3:
                # Sort the atoms into a tuple
                key = (i, j, k)

                # The atoms in forward or reverse order are duplicates
                if key in results or key[::-1] in results:
                    continue

                # Process the filters
                if filter_atoms(*key, **filters):
                    continue

                # Add the distance to the results
                angle = round(measure_angle(i, j, k), 1)
                results[key] = angle

    # Return the sorted list of tuples
    sort_key = lambda i: (i[0][0].chain.id,
                          i[0][1].chain.id,
                          i[0][2].chain.id,
                          i[0][0].residue.number,
                          i[0][1].residue.number,
                          i[0][2].residue.number)
    return [(k[0], k[1], k[2], v)
            for k,v in sorted(results.items(), key=sort_key)]


def measure_dihedrals(molecule, selector1, selector2, selector3, selector4,
                      **filters):
    """Measure the dihedral angles for atoms selected by selector1, selector 2,
    selector3 and selector4.

    Parameters
    ----------
    molecule: :obj:`mollib.Molecule`
        The molecule object from which distances will be measured.
    selector1: str
        A string for the selector of the first atom(s) of the angle.
    selector2: str
        A string for the selector of the second atom(s) of the angle.
    selector3: str
        A string for the selector of the third atom(s) of the angle.
    selector4: str
        A string for the selector of the fourth atom(s) of the angle.
    filters: dict
        See :func:`mollib.core.filter_atoms`


        .. note:: The locator respects the rules outlined in
                  :class:`mollib.Molecule.get_atoms`

    Returns
    -------
    list of tuples
        A sorted list of tuples with (atom1, atom2, atom3, atom4, angle) for
        each angle selected.

    Examples
    --------
    >>> from mollib.core import Molecule, measure_dihedrals
    >>> mol = Molecule('2KXA')
    >>> args = ('3:5.C', '3:5.N', '3:5.CA', '3:5.C')
    >>> kwargs = {'bonded': True}
    >>> dihs = measure_dihedrals(mol, *args, **kwargs)
    >>> for i in dihs: print(i)
    (A.F3.C, A.G4.N, A.G4.CA, A.G4.C, -57.2)
    (A.G4.C, A.A5.N, A.A5.CA, A.A5.C, -68.4)

    """
    # Returned dict
    results = {}

    # This function logs an error if a1 or a2 isn't properly
    # formatted. An additional message is not needed. Just skip
    # it if both atoms aren't found.
    atoms1 = molecule.get_atoms(selector1)
    atoms2 = molecule.get_atoms(selector2)
    atoms3 = molecule.get_atoms(selector3)
    atoms4 = molecule.get_atoms(selector4)

    if not atoms1 or not atoms2 or not atoms3 or not atoms4:
        return []

    for i in atoms1:
        for j in atoms2:
            for k in atoms3:
                for l in atoms4:
                    # Sort the atoms into a tuple
                    key = (i, j, k, l)

                    # The atoms in forward or reverse order are duplicates
                    if key in results or key[::-1] in results:
                        continue

                    # Process the filters
                    if filter_atoms(*key, **filters):
                        continue

                    # Add the distance to the results
                    dihedral = round(measure_dihedral(i, j, k, l), 1)
                    results[key] = dihedral

    # Return the sorted list of tuples
    sort_key = lambda i: (i[0][0].chain.id,
                          i[0][1].chain.id,
                          i[0][2].chain.id,
                          i[0][3].chain.id,
                          i[0][0].residue.number,
                          i[0][1].residue.number,
                          i[0][2].residue.number,
                          i[0][3].residue.number)
    return [(k[0], k[1], k[2], k[3], v)
            for k,v in sorted(results.items(), key=sort_key)]
