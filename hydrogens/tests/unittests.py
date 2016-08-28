import unittest

from mollib.core import Molecule, measure_angle, measure_distance
from mollib.hydrogens import add_hydrogens


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

    def _test_residues(self, molecule_h, molecule_ref, residue_dict, tolerance,
                       skip_first=True):
        """Test that all of the hydrogen atoms listed in residue_dict match
        the same position between molecule_h and molecule_ref, within tolerance.

        Parameters
        ----------
        molecule_h: :obj:`molecule`
            The molecule that has had hydrogens added to it.
        molecule_ref: :obj:`molecule`
            The reference molecule with the correct hydrogen positions
        residue_dict: dict
            A dictionary of residue names (keys) and list of tuples with the
            heavy atom and proton atom names to test.
            ex: {'ARG'] = [('NE', 'HE'), ]}
        tolerance: float
            The tolerance to use for the hydrogen positions between molecule_h
            and molecule_ref to be considered a match.
        """
        for residue in molecule_h.residues:
            if residue.name not in residue_dict:
                continue

            for target_name, h_name in residue_dict[residue.name]:
                # See how far these are from the reference structure
                a1 = residue.get(h_name, None)
                residue_ref = molecule_ref[residue.chain.id][residue.number]
                a2 = residue_ref.get(h_name, None)

                if a1 is None:
                    msg = "molecule_h ({}): ".format(molecule_h.name)
                    msg += "Atom '{}' in residue '{}' ".format(h_name,
                                                               residue_ref)
                    msg += "not found in the hydrogenated molecule."
                    print(msg)
                    continue
                if a2 is None:
                    msg = "molecule_ref ({}): ".format(molecule_ref.name)
                    msg += "Atom '{}' in residue '{}' ".format(h_name,
                                                               residue_ref)
                    msg += "not found in the reference molecule."
                    print(msg)
                    continue
                distance = measure_distance(a1, a2)  # Angstroms

                msg = "Molecule ({}):".format(a1.molecule.name)
                msg += "{} and {} are {:.3f} A ".format(a1, a2, distance)
                msg += "from each other."
                self.assertTrue(within(distance, 0.0, self.tolerance),
                                msg=msg)

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
            mol_ref = Molecule(mol_name)

            mol = Molecule(mol_name)
            add_hydrogens(mol, strip=True)

            self._test_residues(mol, mol_ref, res, self.tolerance)

    def test_add_two_sp2_h(self):
        # These are the sp2 atoms that need two Hs in proteins.
        res = {}
        res['ARG'] = [('NH1', 'HH11'), ('NH1', 'HH12'),
                      ('NH2', 'HH21'), ('NH2', 'HH22')]
        res['ASN'] = [('ND2', 'HD21'), ('ND2', 'HD22')]
        res['GLN'] = [('NE2', 'HE21'), ('NE2', 'HE22')]

        # Reference structure with protons. Try multiple PDBs
        for mol_name in ('2MJB',):  # has RNQ
            mol_ref = Molecule(mol_name)
            mol = Molecule(mol_name)
            add_hydrogens(mol, strip=True)
            self._test_residues(mol, mol_ref, res, self.tolerance)

    def test_add_two_sp3_h(self):
        # These are the sp2 atoms that need two Hs in proteins.
        res = {}
        res['PRO'] = [('CB', 'HB2'), ('CB', 'HB3'),
                      ('CG', 'HG2'), ('CG', 'HG3'),
                      ('CD', 'HD2'), ('CB', 'HD3'), ]
        res['GLY'] = [('CA', 'HA2'), ('CA', 'HA2'), ]
        res['ARG'] = [('CB', 'HB2'), ('CB', 'HB3'),
                      ('CG', 'HG2'), ('CG', 'HG3'),
                      ('CD', 'HD2'), ('CD', 'HD3'), ]
        res['ASN'] = [('CB', 'HB2'), ('CB', 'HB3'), ]
        res['ASP'] = [('CB', 'HB2'), ('CB', 'HB3'), ]
        res['CYS'] = [('CB', 'HB2'), ('CB', 'HB3'), ]
        res['GLN'] = [('CB', 'HB2'), ('CB', 'HB3'),
                      ('CG', 'HG2'), ('CG', 'HG3'), ]
        res['GLU'] = [('CB', 'HB2'), ('CB', 'HB3'),
                      ('CG', 'HG2'), ('CG', 'HG3'), ]
        res['HIS'] = [('CB', 'HB2'), ('CB', 'HB3'), ]
        res['ILE'] = [('CG1', 'HG12'), ('CG1', 'HG13'), ]
        res['LEU'] = [('CB', 'HB2'), ('CB', 'HB3'), ]
        res['LYS'] = [('CB', 'HB2'), ('CB', 'HB3'),
                      ('CG', 'HG2'), ('CG', 'HG3'),
                      ('CD', 'HD2'), ('CD', 'HD3'),
                      ('CE', 'HE2'), ('CE', 'HE3'), ]
        res['MET'] = [('CB', 'HB2'), ('CB', 'HB3'),
                      ('CG', 'HG2'), ('CG', 'HG3'), ]
        res['PHE'] = [('CB', 'HB2'), ('CB', 'HB3'), ]
        res['SER'] = [('CB', 'HB2'), ('CB', 'HB3'), ]
        res['TRP'] = [('CB', 'HB2'), ('CB', 'HB3'), ]
        res['TYR'] = [('CB', 'HB2'), ('CB', 'HB3'), ]

        # Reference structure with protons. Try multiple PDBs
        for mol_name in ('2MJB',):  # has PGRNDQEHILMFSY
                                    # has CW
            mol_ref = Molecule(mol_name)
            mol = Molecule(mol_name)
            add_hydrogens(mol, strip=True)
            self._test_residues(mol, mol_ref, res, self.tolerance)