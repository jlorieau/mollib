cdef class Box:
    """A Box class for storing points in nodes organized in 3-d space and
    quickly returning the nearest neighbor points to a point within a cutoff.
    """

    cdef:
        double x_min, x_max, y_min, y_max, z_min, z_max, _node_size
        int _x_no_box, _y_no_box, _z_no_box
        list _box
        bint preserve_cache_wb_rotation, preserve_cache_wb_translation
        bint preserve_cache_renumber_atoms

    cdef add_points(self, list points)

    cdef get_nodes(self, object point, double radius)