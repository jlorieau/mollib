import unittest

from mollib.core import Molecule, measure_angle, measure_distance
from mollib.hydrogens import add_one_sp2_h, settings


aminoacids = {'ALA', 'GLY', 'SER', 'THR', 'MET', 'CYS', 'ILE', 'LEU',
              'VAL', 'PHE', 'TYR', 'TRP', 'ASN', 'GLN', 'ASP', 'GLU',
              'HIS', 'PRO', 'ARG', 'LYS'}


def within(parameter1, parameter2, tolerance):
    """Tests whether parameter1 and parameter2 are within the tolerance."""
    return abs(parameter1 - parameter2) < tolerance

class TestHydrogenate(unittest.TestCase):
    "Test the hydrogenate functions."

    tolerance = 0.12
    """Tolerance (in Angstroms) to be considered a match to the reference
    structure.


    .. note: The reference structures have been refined against RDCS and the
             H positions may be tilted away from their ideal values. It is for
             this reason that the tolerance is so large.
    """

    def test_add_one_sp2_h(self):
        # TODO: This test only tests atoms that need Hs and have 2 heavy atoms.

        # These are the sp2 atoms that need one H in proteins. These are all
        # sp2 centers with 2 heavy atoms bonded to them.
        res = {}
        res['PHE'] = [('CD1', 'HD1'), ('CD2', 'HD2'), ('CE1', 'HE1'),
                      ('CE2', 'HE2'), ('CZ', 'HZ')]
        res['ARG'] = [('NE', 'HE'), ]
        res['HIS'] = [('ND1', 'HD1'), ('CE1', 'HE1'), ('NE2', 'HE2'),
                      ('CD2', 'HD2')]
        res['TRP'] = [('CD1', 'HD1'), ('NE1', 'HE1'), ('CZ2', 'HZ2'),
                      ('CE3', 'HE3'), ('CZ3', 'HZ3'), ('CH2', 'HH2')]
        res['TYR'] = [('CD1', 'HD1'), ('CD2', 'HD2'), ('CE1', 'HE1'),
                      ('CE2', 'HE2'),]

        # Add amide protons
        for aa in aminoacids:
            if aa == 'PRO':
                continue
            l = res.setdefault(aa, [])
            l.append(('N', 'H'))

        # Reference structure with protons. Try multiple PDBs
        for mol_name in  ('2KXA',  # has FWY
                          '2MJB',  # has FRHY
                          ):
            mol_h = Molecule(mol_name)
            mol = Molecule(mol_name)
            mol.strip_atoms(element='H')

            # Try a variety of atoms. Check that these are within self.tolerance
            # of the reference structure
            for residue in mol.residues:
                if residue.name not in res:  # only check residues from the res
                    continue                 # dict
                for target_name, h_name in res[residue.name]:
                    if residue.first and target_name == 'N':  # NH3+ isn't sp2
                        continue
                    bond_length = (settings.bond_length['N-H']
                                   if target_name.startswith('N') else
                                   settings.bond_length['C-H'])
                    msg = "Could not add {} to {} for {}.".format(h_name,
                                                                  target_name,
                                                                  residue)
                    self.assertTrue(add_one_sp2_h(residue[target_name],
                                                  bond_length), msg=msg)
                    # See how far these are from the reference structure
                    a1 = residue.get(h_name, None)
                    a2 = mol_h['A'][residue.number].get(h_name, None)
                    if a1 is None or a2 is None:
                        continue
                    distance = measure_distance(a1, a2)  # Angstroms
                    print(distance, a1, a2)
                    msg = "Molecule ({}):".format(a1.molecule.name)
                    msg += "{} and {} are {:.3f} A ".format(a1, a2, distance)
                    msg += "from each other."
                    self.assertTrue(within(distance, 0.0, self.tolerance),
                                    msg=msg)
