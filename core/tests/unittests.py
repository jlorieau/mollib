"""Unittests for the core module.
"""
# Author: Justin L Lorieau
# Copyright 2016

import unittest

aminoacids = {'ALA', 'GLY', 'SER', 'THR', 'MET', 'CYS', 'ILE', 'LEU',
              'VAL', 'PHE', 'TYR', 'TRP', 'ASN', 'GLN', 'ASP', 'GLU',
              'HIS', 'PRO', 'ARG', 'LYS'}


from mollib.core import Molecule

class TestMolLib(unittest.TestCase):

    def test_large_molecule(self):
        "Tests the parsing and performance of a very large protein complex."
        import string

        mol = Molecule('3H0G')

        # Test that all of the chains were read in correctly
        chains = list(string.ascii_uppercase)[:24]  # chains A-X
        chains += ['A*', 'B*', 'C*', 'I*', 'J*', 'L*', 'M*', 'N*', 'O*', 'U*',
                   'V*', 'X*']
        self.assertEqual([c.id for c in mol.chains], sorted(chains))
        self.assertEqual(mol.chain_size, len(chains))

        # Test the molecular mass of each chain
        for chain in mol.chains:
            self.assertGreater(chain.mass, 0.)

        self.assertAlmostEqual(mol.mass, 833388.28, 2)

    def test_multiple_models(self):
        """Tests reading PDB files with multiple models. Only the first model
        should be read in."""
        mol = Molecule('2KXA')  # 20 models

        # These are the coordinates for this atom of the first model
        self.assertEqual(mol['A'][3]['N'].pos[0], 13.766)
        self.assertEqual(mol['A'][3]['N'].pos[1], -3.965)
        self.assertEqual(mol['A'][3]['N'].pos[2], 5.893)

    def test_residue_ordering(self):
        """Tests the linked lists of residues."""
        # Molecule is influenza M2 (19-49). It has 4 chains (A, B, C, D) and
        # one drug 11.
        mol = Molecule('2MUV')

        chain_ids = [c.id for c in mol.chains]
        self.assertEqual(chain_ids,
                         ['A', 'B', 'C', 'C*', 'D'])

        for chain_id in chain_ids:
            residue_size = mol[chain_id].residue_size

            # Only the first residue is residue.first
            first = [r.first for r in mol[chain_id].residues]
            self.assertEqual(first,
                             [True] + [False]*(residue_size-1))

            # Only the last residue is residue.last
            last = [r.last for r in mol[chain_id].residues]
            self.assertEqual(last,
                             [False] * (residue_size - 1) + [True])

            # Check the flags and settings for the first and last residue
            for count, residue in enumerate(mol[chain_id].residues):
                if count == 0: # This is the first residue
                    self.assertTrue(residue.first)
                    self.assertIsNone(residue.prev_residue)
                elif count + 1 == residue_size: # This is the last residue
                    self.assertTrue(residue.last)
                    self.assertIsNone(residue.next_residue)

            # Checking the linking. These tests have to use __repr__ because
            # the prev_residue and next_residue are only weakref proxies.
            prev_residues = [r.prev_residue if r is not None else None
                             for r in mol[chain_id].residues]
            self.assertEqual([r.__repr__() for r in prev_residues],
                             ['None'] + [r.__repr__()
                                       for r in mol[chain_id].residues][:-1])

            next_residues = [r.next_residue if r is not None else None
                             for r in mol[chain_id].residues]
            self.assertEqual([r.__repr__() for r in next_residues],
                             [r.__repr__()
                              for r in mol[chain_id].residues][1:] + ['None'])

    def test_protein_topologies(self):
        """Tests whether the atom topologies have been correctly set for
        proteins."""

        # Function to check all of the each atom topologies
        def check_topology(molecule):
            for residue in mol.residues:
                if residue.name not in aminoacids:
                    continue

                # Check the 'N' atom
                if residue.first and residue.name != 'PRO':
                    self.assertEqual(residue['N'].topology,
                                     {'CA', 'H1', 'H2', 'H3'})
                elif residue.name == 'PRO':
                    self.assertEqual(residue['N'].topology,
                                     {'CA', 'CD', 'C-1'})
                else:
                    self.assertEqual(residue['N'].topology,
                                     {'CA', 'H', 'C-1'})

                # Check the 'C', 'O'  and 'OXT atom
                if residue.last:
                    self.assertEqual(residue['C'].topology,
                                     {'CA', 'O', 'OXT'})
                    self.assertIn('C', residue['O'].topology)
                    self.assertEqual(residue['OXT'].topology,
                                     {'C', 'HXT'})
                else:
                    self.assertEqual(residue['C'].topology,
                                     {'CA', 'O', 'N+1'})
                    self.assertIn('C', residue['O'].topology)

                # check the 'CA atom
                if residue.name == 'GLY':
                    self.assertEqual(residue['CA'].topology,
                                     {'N', 'HA2', 'HA3', 'C'})
                else:
                    self.assertEqual(residue['CA'].topology,
                                     {'N', 'CB', 'HA', 'C'})

                # Check the CB
                if residue.name in ('PRO', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN',
                                    'GLU', 'HIS', 'LEU', 'LYS', 'MET', 'PHE',
                                    'SER', 'TRP', 'TYR'):
                    t = residue['CB'].topology
                    test = any((t == {'CA', 'HB2', 'HB3','CG'},
                                t == {'CA', 'HB2', 'HB3', 'SG'},
                                t == {'CA', 'HB2', 'HB3', 'OG'}))
                    self.assertTrue(test)
                elif residue.name is 'ILE':
                    self.assertEqual(residue['CB'].topology,
                                     {'HB', 'CA', 'CG1', 'CG2'})
                elif residue.name is 'THR':
                    self.assertEqual(residue['CB'].topology,
                                     {'HB', 'CA', 'CG2', 'OG1'})
                elif residue.name is 'VAL':
                    self.assertEqual(residue['CB'].topology,
                                     {'CA', 'HB', 'CG1', 'CG2'})

        # Molecule is a domain of trypsinogen with 3 cysteine bridges
        # and 2 Calcium ions
        mol = Molecule('2PTN')

        # There should be 20 CONECT entries
        self.assertEqual(len(mol.connections), 20)

        # The following tests the topology of the cysteine bridges
        pairs = ((22, 157), (42, 58), (128, 232), (136, 201), (168, 182),
                 (191, 220))
        for i,j in pairs:
            # Test the atom topologies
            self.assertEqual(mol['A'][i]['SG'].topology,
                             {'2PTN.A.C{}-SG'.format(j), 'CB'})
            self.assertEqual(mol['A'][j]['SG'].topology,
                             {'2PTN.A.C{}-SG'.format(i), 'CB'})

        check_topology(mol)

        # Check a structure with mulitple chains. The M2 channel has 4 chains
        mol = Molecule('2MUV')
        check_topology(mol)

    def test_add_remove_atoms(self):
        """Tests the add and remove atom functions."""
        mol = Molecule('2KXA')

        # Replace H with a methyl to a tyrosine
        Y22 = mol['A'][22]
        h = Y22['HH']
        self.assertIn(h, Y22['OH'].bonded_atoms())
        mol.del_atom(h)
        self.assertNotIn(h, Y22['OH'].bonded_atoms())

        mol.add_atom(name='CH', pos=(0,0,0), charge=0, element='C',
                     residue=Y22, bonded_atoms=[Y22['OH']],)
        self.assertIn(Y22['CH'], Y22['OH'].bonded_atoms())
        for i in '123':
            mol.add_atom(name='HH'+i, pos=(0, 0, 0), charge=0, element='H',
                         residue=Y22, bonded_atoms=[Y22['CH']], )
        self.assertEqual(Y22['CH'].bonded_atoms(),
                         {Y22['OH'], Y22['HH1'], Y22['HH2'], Y22['HH3']})

        # Replace methyl with an H on the tyrosine
        for name in ('OH', 'HH1', 'HH2', 'HH3'):
            mol.del_atom(Y22[name])



