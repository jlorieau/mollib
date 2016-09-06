import unittest

from mollib.core import Molecule


class TestAtom(unittest.TestCase):

    def test_in_topology(self):
        "Tests that the 'in_topology' function returns the correct value."

        mol = Molecule('2KXA')

        # Test for false topologies
        self.assertFalse(mol['A'][3]['N'].in_topology(mol['A'][4]['H']))

        # Test for actual backbone connectivities
        for res in mol.residues:
            res_num_i = res.number
            chain_i = res.chain.id

            res_num_h = res.number - 1
            res_num_j = res.number + 1

            # Check the previous residue, if it exists
            if (chain_i in mol and res_num_h in mol[chain_i]):
                prev_res = mol[chain_i][res_num_h]
                if 'H' in res:
                    self.assertTrue(all(map(res['N'].in_topology,
                                        [prev_res['C'], res['H'], res['CA']])))
                    self.assertTrue(all(map(res['H'].in_topology,
                                            [res['N'], ])))
                else:
                    self.assertTrue(all(map(res['N'].in_topology,
                                            [prev_res['C'], res['H'],
                                             res['CA']])))
                    self.assertTrue(all(map(res['H'].in_topology,
                                            [res['N'], ])))

            # Check the next residue, if it exists
            if (chain_i in mol and res_num_j in mol[chain_i]):
                next_res = mol[chain_i][res_num_j]
                self.assertTrue(all(map(res['C'].in_topology,
                                        [res['CA'], next_res['N']])))

            self.assertTrue(all(map(res['CA'].in_topology,
                                    [res['N'], res['C']])))

