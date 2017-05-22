"""
The Hydrogen Bonds (hbonds) module is use to identify and characterize hydrogen 
bonds in biomolecules and to classify secondary structure units. 
"""

from .hbonds import find_hbond_partners, find_dipoles
from .classify_secondary_structure import classify_residues
from . import settings

# Load plugins
from .plugin import Hbonds
hbonds = Hbonds()

from mollib.core import register_settings
register_settings(settings)

# TODO: add the '-o' option for printing reports.
