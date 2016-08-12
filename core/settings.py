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

default_pH = 7.0
"""The default pH of new molecules."""

pKs = {'ASP': 3.5,
       'GLU': 4.2,
       'HIS': 6.6,
       'CYS': 6.8,
       'TYR': 10.3,
       'LYS': 10.5,
       'C-term': 3.3,
       'N-term': 7.7}
"""The default pKs of ionizable amino-acids. [Ref]_


  .. [Ref]: G. R. Grimsley, J. M. Scholtz, C. N. Pace,
            Protein Sci. 18, 247-51 (2009).
"""
