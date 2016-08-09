"""
MolLib default settings
"""
# Author: Justin L Lorieau
# Copyright: 2016

# TODO: move plugin settings to their own plugin submodule

# Ramachandran parameters

dihedral_helix_angles = (-60, -45)
"The default optimal Ramachandran angles for a helix (in degrees)."

dihedral_helix_angles_threshold = (30, 30)
"""The default Ramanchandran angle threshold (+/- threshold) to be considered a
helix (in degrees)"""
