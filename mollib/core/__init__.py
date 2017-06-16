"""
The mollib is a simple Python toolset for manipulating molecular
structures in Python. 
"""

from . import settings
from .settings_manager import (register_settings, list_global_settings,
                               load_settings, import_settings)
from . import plugins
from .atom import Atom, sorted_atom_list
from .chain import Chain
from .geometry import *
from .molecule import Molecule
from .residue import Residue
from .readers import MoleculeReader
from .utils import *

# TODO: add a '--model' option to select models other than the first.
