# cython: profile=False

cdef extern from "math.h":
    double sqrt(double m) nogil
    double ceil(double m) nogil

from itertools import chain

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

    # Get an iterable of atoms to search
    if atom_selection is None:
        element_list = elements.split('|') if elements != '' else []

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


cdef class Box:
    """A Box class for storing points in nodes organized in 3-d space and
    quickly returning the nearest neighbor points to a point within a cutoff.
    """

    cdef double x_min, x_max, y_min, y_max, z_min, z_max, _node_size
    cdef int _x_no_box, _y_no_box, _z_no_box
    cdef list _box

    @cython.cdivision(True)
    def __init__(self, point_list, double node_size=2.0):
        assert node_size > 0.0

        cdef double x_min, x_max, y_min, y_max, z_min, z_max
        cdef int i, j, k, no_points
        cdef object point
        cdef list points

        points = list(point_list)  # Convert generator to a list to traverse
                                   # multiple times

        # Get the first point to assign the min/max values
        point = points[0]

        x_min = point[0]
        x_max = point[0]
        y_min = point[1]
        y_max = point[1]
        z_min = point[2]
        z_max = point[2]

        no_points = len(points)
        for i in range(1, no_points):
            point = points[i]

            if x_min > point[0]:
                x_min = point[0]
            if x_max < point[0]:
                x_max = point[0]

            if y_min > point[1]:
                y_min = point[1]
            if y_max < point[1]:
                y_max = point[1]

            if z_min > point[2]:
                z_min = point[2]
            if z_max < point[2]:
                z_max = point[2]

        x_min -= node_size
        y_min -= node_size
        z_min -= node_size

        x_max += node_size
        y_max += node_size
        z_max += node_size

        self.x_min = x_min
        self.x_max = x_max
        self.y_min = y_min
        self.y_max = y_max
        self.z_min = z_min
        self.z_max = z_max

        self._x_no_box = (int) ((x_max - x_min)/node_size)
        self._y_no_box = (int) ((y_max - y_min)/node_size)
        self._z_no_box = (int) ((z_max - z_min)/node_size)

        self._node_size = node_size

        # Initialize the box
        box = list()

        for i in range(self._x_no_box):
            box.append(list())
            for j in range(self._y_no_box):
                box[i].append(list())
                for k in range(self._z_no_box):
                    box[i][j].append(None)

        self._box = box

        # Add points to the box
        self.add_points(points)

    def dist(self, p1, p2):
        "Calculate the cutoff between points p1 and p2"
        cdef float dist = 0.0
        cdef int i
        cdef point_size = len(p1)
        for i in range(point_size):
            dist += (p1[i] - p2[i])**2
        return sqrt(dist)

    def size(self):
        return len(self._box)*len(self._box[0])*len(self._box[0][0])

    def num_points(self):
        count = 0
        for i in range(len(self._box)):
            for j in range(len(self._box[i])):
                    count += sum([len(k) for k in self._box[i][j] if k is not None])
        return count

    cdef add_points(self, list points):
        "Add points to the box."
        cdef int i = 0
        cdef int no_points = len(points)

        for i in range(no_points):
            point = points[i]
            node_set = self.get_node(point, create=True)
            node_set.append(point)

    def get_node(self, point, create=False):
        "Return the box node at the given point."
        cdef int int_x, int_y, int_z
        int_x = (int) ((point[0] - self.x_min) / self._node_size)
        int_y = (int) ((point[1] - self.y_min) / self._node_size)
        int_z = (int) ((point[2] - self.z_min) / self._node_size)

        if create and self._box[int_x][int_y][int_z] is None:
            self._box[int_x][int_y][int_z] = list()

        return self._box[int_x][int_y][int_z]

    cdef get_nodes(self, object point, double radius):
        cdef int x_min, x_max, y_min, y_max, z_min, z_max, i, j
        cdef list results

        x_min = (int) ((point[0] - radius - self.x_min) / self._node_size)
        x_max = (int) (ceil((point[0] + radius - self.x_min) / self._node_size))

        x_max = x_min + 1 if x_min == x_max else x_max
        x_max = max(x_max, 0)
        x_max = min(x_max, self._x_no_box)
        x_min = min(x_min, self._x_no_box)
        x_min = max(x_min, 0)

        y_min = (int) ((point[1] - radius - self.y_min) / self._node_size)
        y_max = (int) (ceil((point[1] + radius - self.y_min) / self._node_size))

        y_max = y_min + 1 if y_min == y_max else y_max
        y_max = max(y_max, 0)
        y_max = min(y_max, self._y_no_box)
        y_min = min(y_min, self._y_no_box)
        y_min = max(y_min, 0)

        z_min = (int) ((point[2] - radius - self.z_min) / self._node_size)
        z_max = (int) (ceil((point[2] + radius - self.z_min) / self._node_size))

        z_max = z_min + 1 if z_min == z_max else z_max
        z_max = max(z_max, 0)
        z_max = min(z_max, self._z_no_box)
        z_min = min(z_min, self._z_no_box)
        z_min = max(z_min, 0)

        results = []
        for i in range(x_min, x_max):
            for j in range(y_min, y_max):
                results += self._box[i][j][z_min:z_max]

        return filter(None, results)

    def get_points(self, point, radius):
        nodes = self.get_nodes(point, radius)

        return (a for a in chain(*nodes) if self.dist(a, point) < radius)
