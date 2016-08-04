"""
Utility functions.

   @Author:             Justin L Lorieau <jlorieau>
   @Date:               2016-08-04T12:07:09-05:00
   @Last modified by:   jlorieau
   @Last modified time: 2016-08-04T12:42:04-05:00
   @License:            Copyright 2016
"""
from math import sqrt


def vector_length(vector):
    "Returns the length (in A) of a vector"
    return sqrt(sum([i*i for i in vector]))


def calc_vector(atom_i, atom_j, normalize=True):
    "Returns the vector between atoms 'i' and 'j' with optional normalization."
    vec = atom_i.pos - atom_j.pos

    if normalize:
        length = vector_length(vec)
        return vec / length
    else:
        return vec
