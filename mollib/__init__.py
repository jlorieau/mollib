from .core import Atom, Residue, Chain, Molecule
from . import hydrogens
from . import hbonds
from . import statistics
from . import pa
from .__main__ import main

from mollib.utils.version import get_version

VERSION = (1, 0, 0, 'alpha', 1)

__version__ = get_version(VERSION)
