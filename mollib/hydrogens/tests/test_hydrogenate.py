import unittest

from mollib.core import (Molecule, measure_angle, measure_distance,
                         measure_dihedral)
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
        for mol_name in ('2MJB',   # has PGRNDQEHILMFSY
                         '2A5M'):  # has CW
            mol_ref = Molecule(mol_name)
            mol = Molecule(mol_name)
            add_hydrogens(mol, strip=True)
            self._test_residues(mol, mol_ref, res, self.tolerance)

    def test_add_three_sp3_h(self):
        # These are the sp2 atoms that need two Hs in proteins.
        res = {}
        res['ALA'] = [('CB', 'HB1'), ('CB', 'HB2'), ('CB', 'HB3')]
        res['ILE'] = [('CG2', 'HG21'), ('CG2', 'HG22'), ('CG2', 'HG23'),
                      ('CD1', 'HD11'), ('CD1', 'HD12'), ('CD1', 'HD13')]
        res['LEU'] = [('CD1', 'HD11'), ('CD1', 'HD12'), ('CD1', 'HD13'),
                      ('CD2', 'HD21'), ('CD2', 'HD22'), ('CD2', 'HD23')]
        res['LYS'] = [('NZ', 'HZ1'), ('NZ', 'HZ2'), ('NZ', 'HZ3')]
        res['MET'] = [('CE', 'HE1'), ('CE', 'HE2'), ('CE', 'HE3')]
        res['THR'] = [('CG2', 'HG21'), ('CG2', 'HG22'), ('CG2', 'HG23')]
        res['VAL'] = [('CG1', 'HG11'), ('CG1', 'HG12'), ('CG1', 'HG13'),
                      ('CG2', 'HG21'), ('CG2', 'HG22'), ('CG2', 'HG23')]


        angle_offsets = {('ALA', 'CB'): 180.,
                         ('ILE', 'CD1'): 180.,
                         ('ILE', 'CG2'): 180.,
                         ('LEU', 'CD1'): 180.,
                         ('LEU', 'CD2'): 180.,
                         ('LYS', 'NZ'): 180.,
                         ('MET', 'CE'): 180.,
                         ('THR', 'CG2'): 180.,
                         ('VAL', 'CG1'): 180.,
                         ('VAL', 'CG2'): 180.,}
        # Reference structure with protons. Try multiple PDBs
        for mol_name in ('2MJB',):  # has AILKMTV
            mol_ref = Molecule(mol_name)
            mol = Molecule(mol_name)

            # First set the alpha angle for all relevant sp3 atoms so that
            # the atom positions match with methyl or amine rotations.
            for residue in mol_ref.residues:
                if residue.name in res:
                    for atom_name, h_name in res[residue.name]:
                        if not h_name.endswith('1'):
                            continue

                        atom, h = residue[atom_name], residue[h_name]
                        bonded = atom.bonded_heavy_atoms(sorted=True)[0]
                        bonded2 = [a for a
                                   in bonded.bonded_heavy_atoms(sorted=True)
                                   if a != atom][0]

                        # The angles in the PDB files are different from these
                        # dihedrals by angle offsets.
                        offset = angle_offsets.get((residue.name, atom_name),
                                                   0.0)
                        angle = (measure_dihedral(h, atom, bonded, bonded2)
                                 + offset)

                        mol.set_parameter('Add_hydrogens',
                                          atom.fullname + '_alpha',
                                          angle)

            add_hydrogens(mol, strip=True)
            self._test_residues(mol, mol_ref, res, self.tolerance)

    def test_ionizeable_groups(self):
        """Tests the correct protonation of ionizeable groups.
        """

        def test_geometry(atom, expected_no_H, expected_angle=None,
                          tolerance=6):
            "Tests an atom for the expect number of hydrogens and geometry."
            # Count the number of hydrogens
            hydrogens = [b for b in atom.bonded_atoms() if b.element == 'H']

            msg = "Atom '{}' has {} hydrogens. {} expected."
            msg = msg.format(atom, len(hydrogens), expected_no_H)
            self.assertEqual(len(hydrogens), expected_no_H, msg)

            # Check the angles and geometry
            if expected_no_H > 0:
                bonded_heavy_atom = atom.bonded_heavy_atoms()[0]
                angle = measure_angle(bonded_heavy_atom, atom, hydrogens[0])

                msg = "The '{}'-'{}'-'{}' angle is {}. Expected {}."
                msg = msg.format(bonded_heavy_atom, atom, hydrogens[0],
                                 angle, expected_angle)
                self.assertTrue(within(angle, expected_angle, tolerance),
                                msg)

        def test_either_geometry(atom_1, atom_2, expected_no_H,
                                 expected_angle=None, tolerance=6):
            """Tests whether atom_1 and/or atom_2 has the expected number of
            hydrogens. However, atom_1 and atom_2 can only have one hydrogen
            at a time. Also checks geometry
            """
            hydrogens_1 = [b for b in atom_1.bonded_atoms()
                           if b.element == 'H']
            hydrogens_2 = [b for b in atom_2.bonded_atoms()
                           if b.element == 'H']

            msg = "Atom '{}' has {} hydrogens. <= 1 expected."
            self.assertLessEqual(len(hydrogens_1), 1,
                                 msg=msg.format(atom_1, len(hydrogens_1)))
            self.assertLessEqual(len(hydrogens_2), 1,
                                 msg=msg.format(atom_2, len(hydrogens_2)))

            msg = "Atoms '{}' and '{}' have {} hydrogens together. {} expected."
            total_hydrogens = len(hydrogens_1)+len(hydrogens_2)
            self.assertEqual(total_hydrogens,
                             expected_no_H,
                             msg=msg.format(atom_1, atom_2, total_hydrogens,
                                            expected_no_H))

            if expected_no_H > 0:
                msg = "The '{}'-'{}'-'{}' angle is {}. Expected {}."

                if len(hydrogens_1) > 0:
                    bonded_heavy_atom = atom_1.bonded_heavy_atoms()[0]

                    angle = measure_angle(bonded_heavy_atom, atom_1,
                                          hydrogens_1[0])
                    self.assertTrue(within(angle, expected_angle, tolerance),
                                    msg=msg.format(bonded_heavy_atom, atom_1,
                                                   hydrogens_1[0], angle,
                                                   expected_angle))
                if len(hydrogens_2) > 0:
                    bonded_heavy_atom = atom_2.bonded_heavy_atoms()[0]

                    angle = measure_angle(bonded_heavy_atom, atom_2,
                                          hydrogens_2[0])
                    self.assertTrue(within(angle, expected_angle, tolerance),
                                    msg=msg.format(bonded_heavy_atom, atom_2,
                                                   hydrogens_2[0], angle,
                                                   expected_angle))

        # Ionizeable groups in proteins: DEHCYK, first residue, last
        mol = Molecule('2MJB')  # has DEHYK

        # test first at high pH
        mol.pH = 12
        add_hydrogens(mol, strip=True)

        for residue in mol.residues:
            if residue.first:
                test_geometry(atom=residue['N'], expected_no_H=2,
                              expected_angle=120.)
            if residue.name == 'ASP':
                test_either_geometry(atom_1=residue['OD1'],
                                     atom_2=residue['OD2'],
                                     expected_no_H=0)
            if residue.name == 'GLU':
                test_either_geometry(atom_1=residue['OE1'],
                                     atom_2=residue['OE2'],
                                     expected_no_H=0)
            if residue.name == 'HIS':
                test_either_geometry(atom_1=residue['ND1'],
                                     atom_2=residue['NE2'],
                                     expected_no_H=1,
                                     expected_angle=120.)
            if residue.name == 'TYR':
                test_geometry(atom=residue['OH'], expected_no_H=0)
            if residue.name == 'LYS':
                test_geometry(atom=residue['NZ'], expected_no_H=2,
                              expected_angle=120.)
            if residue.last:
                test_either_geometry(atom_1=residue['O'],
                                     atom_2=residue['OXT'],
                                     expected_no_H=0)

        # test first at intermediate pH
        mol.pH = 6
        add_hydrogens(mol, strip=True)

        for residue in mol.residues:
            if residue.first:
                test_geometry(atom=residue['N'], expected_no_H=3,
                              expected_angle=109.5)
            if residue.name == 'ASP':
                test_either_geometry(atom_1=residue['OD1'],
                                     atom_2=residue['OD2'],
                                     expected_no_H=0)
            if residue.name == 'GLU':
                test_either_geometry(atom_1=residue['OE1'],
                                     atom_2=residue['OE2'],
                                     expected_no_H=0)
            if residue.name == 'HIS':
                test_either_geometry(atom_1=residue['ND1'],
                                     atom_2=residue['NE2'],
                                     expected_no_H=2,
                                     expected_angle=120.)
            if residue.name == 'TYR':
                test_geometry(atom=residue['OH'], expected_no_H=1,
                              expected_angle=109.5)
            if residue.name == 'LYS':
                test_geometry(atom=residue['NZ'], expected_no_H=3,
                              expected_angle=109.5)
            if residue.last:
                test_either_geometry(atom_1=residue['O'],
                                     atom_2=residue['OXT'],
                                     expected_no_H=0)

        # test first at intermediate pH
        mol.pH = 1
        add_hydrogens(mol, strip=True)

        for residue in mol.residues:
            if residue.first:
                test_geometry(atom=residue['N'], expected_no_H=3,
                              expected_angle=109.5)
            if residue.name == 'ASP':
                test_either_geometry(atom_1=residue['OD1'],
                                     atom_2=residue['OD2'],
                                     expected_no_H=1,
                                     expected_angle=109.5)
            if residue.name == 'GLU':
                test_either_geometry(atom_1=residue['OE1'],
                                     atom_2=residue['OE2'],
                                     expected_no_H=1,
                                     expected_angle=109.5)
            if residue.name == 'HIS':
                test_either_geometry(atom_1=residue['ND1'],
                                     atom_2=residue['NE2'],
                                     expected_no_H=2,
                                     expected_angle=120.)
            if residue.name == 'TYR':
                test_geometry(atom=residue['OH'], expected_no_H=1,
                              expected_angle=109.5)
            if residue.name == 'LYS':
                test_geometry(atom=residue['NZ'], expected_no_H=3,
                              expected_angle=109.5)
            if residue.last:
                test_either_geometry(atom_1=residue['O'],
                                     atom_2=residue['OXT'],
                                     expected_no_H=1,
                                     expected_angle=109.5)
