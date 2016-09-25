from .hbonds import find_hbond_partners, find_dipoles
from .classify_secondary_structure import classify_residues
from . import settings

# Load plugins
from .plugin import Hbonds
hbonds = Hbonds()

from mollib.core import register_settings
register_settings(settings)