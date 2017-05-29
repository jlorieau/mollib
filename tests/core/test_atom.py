import unittest

from mollib.core import Molecule
from mollib.core.atom import re_atom


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

    def test_atom_name_regex(self):
        """Test the matching of the re_atom regex."""

        def match_result(label, expected_result):
            match = re_atom.match(label)
            if expected_result is None:
                msg = ("Did not expect the label '{}' to match, but got the "
                       "following dict instead: {}.")
                if hasattr(match, 'groupdict'):
                    msg = msg.format(label, match.groupdict())
                else:
                    msg = msg.format(label, '')
                self.assertIsNone(match, msg=msg)
                return None

            msg = ("Label '{}' did not match the expected result: {}. Instead "
                   "the following dict was received: {}.")
            self.assertEqual(match.groupdict(), expected_result,
                             msg=msg.format(label, expected_result,
                                            match.groupdict()))

        # The following are valid labels
        label = "13CA"
        match_result(label, {'molecule': None, 'chain_id': None,
                             'residue_letter': None, 'residue_number': '13',
                             'atom_name': 'CA'})


        label = "13.HA2"
        match_result(label, {'molecule': None, 'chain_id': None,
                             'residue_letter': None, 'residue_number': '13',
                             'atom_name': 'HA2'})


        label = "I13.HA2"
        match_result(label, {'molecule': None, 'chain_id': None,
                             'residue_letter': 'I', 'residue_number': '13',
                             'atom_name': 'HA2'})

        label = "A.I13.HA2"
        match_result(label, {'molecule': None, 'chain_id': 'A',
                             'residue_letter': 'I', 'residue_number': '13',
                             'atom_name': 'HA2'})

        label = "2KXA:A.I13.HA2"
        match_result(label, {'molecule': '2KXA', 'chain_id': 'A',
                             'residue_letter': 'I', 'residue_number': '13',
                             'atom_name': 'HA2'})

        # The following are invalid labels
        label = "I13-HA2"
        match_result(label, None)

        label = "13.2HA"
        match_result(label, None)

        # Now try to see if the finditer works to match multiple labels.
        text = """
        Atom Names  Values
        ----------  ------
        I13CA       2.324
        W14N        0.123
        T15C        -2.3e-4
        """

        items = [i.group(0).strip() for i in re_atom.finditer(text)]
        self.assertEqual(items, ['I13CA', 'W14N', 'T15C'])




