"""
Readers for Protein Databank Files.
"""

from .base import MoleculeReader


class PDBReader(MoleculeReader):
    pass

class FuzzyPDBReader(MoleculeReader):
    pass