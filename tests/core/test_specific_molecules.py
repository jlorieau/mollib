import unittest

from mollib import Molecule


class TestMolecules(unittest.TestCase):
    "Test the loading of specific molecules"

    def test_1C1D(self):
        mol = Molecule('1C1D')

        # Check that the chains were properly loaded
        self.assertEqual([c.id for c in mol.chains],
                         ['A', 'A*', 'B', 'B*'])

        # Check that the residues were correctly linked
        # Check chain 'A' and 'B'
        self.assertTrue(mol['A'][1].first)
        self.assertFalse(mol['A'][1].last)
        self.assertIsNone(mol['A'][1].prev_residue)
        self.assertEqual(mol['A'][1].next_residue,
                         mol['A'][2])

        self.assertFalse(mol['A'][349].first)
        self.assertTrue(mol['A'][349].last)
        self.assertEqual(mol['A'][349].prev_residue,
                         mol['A'][348])
        self.assertIsNone(mol['A'][349].next_residue)

        self.assertTrue(mol['B'][1].first)
        self.assertFalse(mol['B'][1].last)
        self.assertIsNone(mol['B'][1].prev_residue)
        self.assertEqual(mol['B'][1].next_residue,
                         mol['B'][2])

        self.assertFalse(mol['B'][348].first)
        self.assertTrue(mol['B'][348].last)
        self.assertEqual(mol['B'][348].prev_residue,
                         mol['B'][347])
        self.assertIsNone(mol['B'][348].next_residue)

        # Check the linking of the HETATM chains
        self.assertTrue(mol['A*'][850].first)
        self.assertTrue(mol['A*'][850].last)
        self.assertIsNone(mol['A*'][850].prev_residue)
        self.assertIsNone(mol['A*'][850].next_residue)

        self.assertTrue(mol['A*'][853].first)
        self.assertTrue(mol['A*'][853].last)
        self.assertIsNone(mol['A*'][853].prev_residue)
        self.assertIsNone(mol['A*'][853].next_residue)

        self.assertTrue(mol['A*'][361].first)
        self.assertTrue(mol['A*'][361].last)
        self.assertIsNone(mol['A*'][361].prev_residue)
        self.assertIsNone(mol['A*'][361].next_residue)

        self.assertTrue(mol['A*'][360].first)
        self.assertTrue(mol['A*'][360].last)
        self.assertIsNone(mol['A*'][360].prev_residue)
        self.assertIsNone(mol['A*'][360].next_residue)

        self.assertTrue(mol['B*'][851].first)
        self.assertTrue(mol['B*'][851].last)
        self.assertIsNone(mol['B*'][851].prev_residue)
        self.assertIsNone(mol['B*'][851].next_residue)

        self.assertTrue(mol['B*'][852].first)
        self.assertTrue(mol['B*'][852].last)
        self.assertIsNone(mol['B*'][852].prev_residue)
        self.assertIsNone(mol['B*'][852].next_residue)

        self.assertTrue(mol['B*'][880].first)
        self.assertTrue(mol['B*'][880].last)
        self.assertIsNone(mol['B*'][880].prev_residue)
        self.assertIsNone(mol['B*'][880].next_residue)

        self.assertTrue(mol['B*'][761].first)
        self.assertTrue(mol['B*'][761].last)
        self.assertIsNone(mol['B*'][761].prev_residue)
        self.assertIsNone(mol['B*'][761].next_residue)

        self.assertTrue(mol['B*'][760].first)
        self.assertTrue(mol['B*'][760].last)
        self.assertIsNone(mol['B*'][760].prev_residue)
        self.assertIsNone(mol['B*'][760].next_residue)

        self.assertTrue(mol['B*'][860].first)
        self.assertTrue(mol['B*'][860].last)
        self.assertIsNone(mol['B*'][860].prev_residue)
        self.assertIsNone(mol['B*'][860].next_residue)


        # Check access to bonded atoms
        self.assertEqual(mol['A*'][361]['N'].bonded_atoms(sorted=True),
                         [mol['A*'][361]['CA']])
        self.assertEqual(mol['A*'][361]['C'].bonded_atoms(sorted=True),
                         [mol['A*'][361]['O'], mol['A*'][361]['CA']])
        self.assertEqual(mol['A*'][361]['CA'].bonded_atoms(sorted=True),
                         [mol['A*'][361]['N'], mol['A*'][361]['C'],
                          mol['A*'][361]['CB']])

