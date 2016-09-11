# cython: profile=False

cdef extern from "math.h":
    double sqrt(double m)

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
    cdef i, size
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
