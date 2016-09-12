# cython: profile=False

cdef extern from "math.h":
    double sqrt(double m) nogil

import numpy as np
cimport numpy as np
cimport cython


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef double vector_length(np.ndarray[np.float64_t, ndim=1] vector):
    """Returns the length (in A) of a vector"""
    cdef Py_ssize_t i
    cdef Py_ssize_t vector_size = vector.shape[0]
    cdef double v2 = 0.0
    for i in range(vector_size):
        v2 += vector[i]*vector[i]
    return sqrt(v2)


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
    >>> from mollib import Molecule
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

cpdef list within_distance(object atom, double cutoff, str elements='',
                           bint exclude_intraresidue=False,
                           atom_selection=None):
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
    """
    # prefetching atoms in a list reduces the execution time from 26s to 5s

    cdef double distance
    cdef list atom_list = []
    cdef list element_list
    cdef double [:] v1, v2
    cdef object a

    element_list = elements.split('|') if elements != '' else []

    # Get an iterable of atoms to search
    if atom_selection is None:
        atoms = atom.molecule.atoms
        # Filter the atoms in the iterable
        atoms = [a for a in atoms if
                 ((a != atom) and
                  (not element_list or a.element in element_list) and
                  (not exclude_intraresidue or a.residue != atom.residue))]

    else:
        atoms = atom_selection

    for a in atoms:
        v1 = atom.pos
        v2 = a.pos

        with nogil:
            distance = _within_distance(v1, v2, cutoff)

        if distance > 0.0:
           atom_list.append((a, distance))

    return atom_list
