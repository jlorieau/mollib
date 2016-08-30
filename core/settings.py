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

pKs = {'ASP': {'OD1-OD2': (-1.0, 3.5)},
       'GLU': {'OE1-OE2': (-1.0, 4.2)},
       'HIS': {'ND1-NE2': (6.6, 14.0)},
       'CYS': {'SG': (6.8,)},
       'TYR': {'OH': (10.3,)},
       'LYS': {'NZ': (10.5, 14.0, 14.0)},
       'last': {'O-OXT': (-1.0, 3.3)},
       'first': {'N': (7.7, 14.0, 14.0),}
     }
"""The default pKs of ionizable amino-acids. [Ref]_

Some amino-acids have degenerate ionizeable atoms; these are listed and
separated by '-' characters. The different pKs for each ionization is
listed in the subsequent items in the tuple.


  .. [Ref]: G. R. Grimsley, J. M. Scholtz, C. N. Pace,
            Protein Sci. 18, 247-51 (2009).
"""
