"""Pytests for the core module molecules.
"""


import unittest
import weakref

from mollib.core import Molecule

aminoacids = {'ALA', 'GLY', 'SER', 'THR', 'MET', 'CYS', 'ILE', 'LEU',
              'VAL', 'PHE', 'TYR', 'TRP', 'ASN', 'GLN', 'ASP', 'GLU',
              'HIS', 'PRO', 'ARG', 'LYS'}


def test_get_weakrefs():
    "Tests the get_weakrefs Molecule class method."
    # First the Molecule._instances list has to be cleared
    Molecule._instances = []

    molecules = {}
    for name in ('2KXA', '2MUV'):
        # No molecule instance exists yet
        assert Molecule.get_weakref(name) is None

        # Create the molecule instance. A dict is used so that the molecule
        # objects aren't garbage collected between iterations of this loop.
        molecules[name] = Molecule(name)

        # Now a molecule instance should exist, and this function should
        # return a weakref to it.
        ref = Molecule.get_weakref(name)
        assert ref is not None
        assert isinstance(ref, weakref.ref)

        # Make sure that the weakref is indeed for that molecule.
        assert ref() == molecules[name]

    # Now try making a duplicate
    mol = Molecule('2KXA')
    ref = Molecule.get_weakref('2KXA')
    assert ref() == molecules['2KXA']
    assert ref() == mol


def test_large_molecule():
    """Tests the parsing and performance of a very large protein complex.

    .. note:: This test takes about 3.3s. The cProfile bottlenecks are:
              - 1.093s : Molecule._match_atom
              - 0.626s : Primitives.__init__
              - 0.320s : bzip.readlines
              - 0.225s : collections.item
              - 0.170s : Molecule.prev_residue
              - 0.165s : Molecule.next_residue
              - 0.101s : collections.__iter__
    """
    import string

    mol = Molecule('3H0G')

    # Test that all of the chains were read in correctly
    chains = list(string.ascii_uppercase)[:24]  # chains A-X
    chains += ['A*', 'B*', 'C*', 'I*', 'J*', 'L*', 'M*', 'N*', 'O*', 'U*',
               'V*', 'X*']
    assert [c.id for c in mol.chains] == sorted(chains)
    assert mol.chain_size == len(chains)

    # Test the molecular mass of each chain
    for chain in mol.chains:
        assert chain.mass > 0.

    assert (mol.mass - 833388.28) < 0.1


def test_multiple_models():
    """Tests reading PDB files with multiple models. Only the first model
    should be read in."""
    mol = Molecule('2KXA')  # 20 models

    # These are the coordinates for this atom of the first model
    assert mol['A'][3]['N'].pos[0] == 13.766
    assert mol['A'][3]['N'].pos[1] == -3.965
    assert mol['A'][3]['N'].pos[2] == 5.893


def test_residue_ordering():
    """Tests the linked lists of residues."""
    # Molecule is influenza M2 (19-49). It has 4 chains (A, B, C, D) and
    # one drug 11.
    mol = Molecule('2MUV')

    chain_ids = [c.id for c in mol.chains]
    assert chain_ids == ['A', 'B', 'C', 'C*', 'D']

    for chain_id in chain_ids:
        residue_size = mol[chain_id].residue_size

        # Only the first residue is residue.first
        first = [r.first for r in mol[chain_id].residues]
        assert first == [True] + [False]*(residue_size-1)

        # Only the last residue is residue.last
        last = [r.last for r in mol[chain_id].residues]
        assert last == [False] * (residue_size - 1) + [True]

        # Check the flags and settings for the first and last residue
        for count, residue in enumerate(mol[chain_id].residues):
            if count == 0: # This is the first residue
                assert residue.first is True
                assert residue.prev_residue is None
            elif count + 1 == residue_size: # This is the last residue
                assert residue.last is True
                assert residue.next_residue is None

        # Checking the linking. These tests have to use __repr__ because
        # the prev_residue and next_residue are only weakref proxies.
        prev_residues = [r.prev_residue if r is not None else None
                         for r in mol[chain_id].residues]
        assert ([r.__repr__() for r in prev_residues] ==
                ['None'] + [r.__repr__()
                            for r in mol[chain_id].residues][:-1])

        next_residues = [r.next_residue if r is not None else None
                         for r in mol[chain_id].residues]
        assert ([r.__repr__() for r in next_residues] ==
                [r.__repr__() for r in mol[chain_id].residues][1:] +
                ['None'])


def test_protein_topologies():
    """Tests whether the atom topologies have been correctly set for
    proteins."""

    # Function to check all of the each atom topologies
    def check_topology(molecule):
        for residue in mol.residues:
            if residue.name not in aminoacids:
                continue

            # Check the 'N' atom
            if residue.first and residue.name != 'PRO':
                assert residue['N'].topology == {'CA', 'H1', 'H2', 'H3'}
            elif residue.name == 'PRO':
                assert residue['N'].topology == {'CA', 'CD', 'C-1'}
            else:
                assert residue['N'].topology == {'CA', 'H', 'C-1'}

            # Check the 'C', 'O'  and 'OXT atom
            if residue.last:
                assert residue['C'].topology == {'CA', 'O', 'OXT'}
                assert 'C' in residue['O'].topology
                assert residue['OXT'].topology == {'C', 'HXT'}
            else:
                assert residue['C'].topology == {'CA', 'O', 'N+1'}
                assert 'C' in residue['O'].topology

            # check the 'CA atom
            if residue.name == 'GLY':
                assert residue['CA'].topology == {'N', 'HA2', 'HA3', 'C'}
            else:
                assert residue['CA'].topology == {'N', 'CB', 'HA', 'C'}

            # Check the CB
            if residue.name in ('PRO', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN',
                                'GLU', 'HIS', 'LEU', 'LYS', 'MET', 'PHE',
                                'SER', 'TRP', 'TYR'):
                t = residue['CB'].topology
                test = any((t == {'CA', 'HB2', 'HB3','CG'},
                            t == {'CA', 'HB2', 'HB3', 'SG'},
                            t == {'CA', 'HB2', 'HB3', 'OG'}))
                assert test is True
            elif residue.name is 'ILE':
                assert residue['CB'].topology == {'HB', 'CA', 'CG1', 'CG2'}
            elif residue.name is 'THR':
                assert residue['CB'].topology == {'HB', 'CA', 'CG2', 'OG1'}
            elif residue.name is 'VAL':
                assert residue['CB'].topology == {'CA', 'HB', 'CG1', 'CG2'}

    # Molecule is a domain of trypsinogen with 3 cysteine bridges
    # and 2 Calcium ions
    mol = Molecule('2PTN')

    # There should be 20 CONECT entries
    assert len(mol.connections) == 20

    # The following tests the topology of the cysteine bridges
    pairs = ((22, 157), (42, 58), (128, 232), (136, 201), (168, 182),
             (191, 220))
    for i,j in pairs:
        # Test the atom topologies
        assert mol['A'][i]['SG'].topology == {'2PTN:A.C{}.SG'.format(j),
                                              'CB'}
        assert mol['A'][j]['SG'].topology == {'2PTN:A.C{}.SG'.format(i),
                                              'CB'}

    check_topology(mol)

    # Check a structure with mulitple chains. The M2 channel has 4 chains
    mol = Molecule('2MUV')
    check_topology(mol)


def test_add_remove_atoms():
    """Tests the add and remove atom functions."""
    mol = Molecule('2KXA')

    # Replace H with a methyl to a tyrosine
    Y22 = mol['A'][22]
    h = Y22['HH']
    assert h in Y22['OH'].bonded_atoms()
    mol.del_atom(h)
    assert h not in Y22['OH'].bonded_atoms()

    mol.add_atom(name='CH', pos=(0,0,0), charge=0, element='C',
                 residue=Y22, bonded_atoms=[Y22['OH']],)
    assert Y22['CH'] in Y22['OH'].bonded_atoms()
    for i in '123':
        mol.add_atom(name='HH'+i, pos=(0, 0, 0), charge=0, element='H',
                     residue=Y22, bonded_atoms=[Y22['CH']], )
    assert (Y22['CH'].bonded_atoms(sorted=True) ==
            [Y22['OH'], Y22['HH3'], Y22['HH2'], Y22['HH1']])

    # Replace methyl with an H on the tyrosine
    for name in ('OH', 'HH1', 'HH2', 'HH3'):
        mol.del_atom(Y22[name])


def test_hetatm_connectivities():
    "Tests that the CONECT records are correctly read in HETATMS."
    # Influenza M2 structure with adamantane-like molecule bound.
    mol = Molecule('2MUV')

    # the molecule should have 39 entries
    assert len(mol.connections) == 39

    het = mol['C*'][100]
    assert het['BR'].bonded_atoms(sorted=True) == [het['C1'],]
    assert het['S'].bonded_atoms(sorted=True) == [het['C1'], het['C12']]
    assert het['C1'].bonded_atoms(sorted=True) == [het['BR'], het['S'],
                                                   het['C2']]
    assert het['C2'].bonded_atoms(sorted=True) == [het['C1'], het['C3'],
                                                   het['H2']]
    assert het['C3'].bonded_atoms(sorted=True) == [het['C12'], het['C2'],
                                                   het['H3']]
    assert het['C12'].bonded_atoms(sorted=True) == [het['S'], het['C5'],
                                                    het['C3']]
    assert het['C5'].bonded_atoms(sorted=True) == [het['N2'], het['C12'],
                                                   het['H5A'], het['H5']]
    assert het['N2'].bonded_atoms(sorted=True) == [het['C10'], het['C5'],
                                                   het['HN2A'], het['HN2']]
    assert het['C10'].bonded_atoms(sorted=True) == [het['N2'], het['C63'],
                                                    het['C62'], het['C61']]
    assert het['C61'].bonded_atoms(sorted=True) == [het['C10'], het['C71'],
                                                    het['H61A'], het['H61']]
    assert het['C62'].bonded_atoms(sorted=True) == [het['C10'], het['C72'],
                                                    het['H62A'], het['H62']]
    assert het['C63'].bonded_atoms(sorted=True) == [het['C10'], het['C73'],
                                                    het['H63A'], het['H63']]
    assert het['C71'].bonded_atoms(sorted=True) == [het['C83'], het['C81'],
                                                    het['C61'], het['H71']]
    assert het['C72'].bonded_atoms(sorted=True) == [het['C82'], het['C81'],
                                                    het['C62'], het['H72']]
    assert het['C73'].bonded_atoms(sorted=True) == [het['C83'], het['C82'],
                                                    het['C63'], het['H73']]
    assert het['C81'].bonded_atoms(sorted=True) == [het['C72'], het['C71'],
                                                    het['H81A'], het['H81']]
    assert het['C82'].bonded_atoms(sorted=True) == [het['C73'], het['C72'],
                                                    het['H82A'], het['H82']]
    assert het['C83'].bonded_atoms(sorted=True) == [het['C73'], het['C71'],
                                                    het['H83A'], het['H83']]


def test_cache():
    "Tests the cache_clear method."
    class Mock(object):
        "A Mock cache object."
        pass

    class TesterMethod(object):
        "A tester for methods and cache clearing scopes"

        def __init__(self, **kwargs):
            for k,v in kwargs.items():
                setattr(self, k, v)

    # A molecule object to add caches to
    mol = Molecule('2KXA')

    # A list of tester methods
    testers = [TesterMethod(scope='preserve_cache_wb_rotation',
                            method='rotate_zyz',
                            args=(0., 0., 0.)),
               TesterMethod(scope='preserve_cache_wb_translation',
                            method='center'),
               TesterMethod(scope='preserve_cache_add_atoms',
                            method='add_atom',
                            args=('T', (0., 0., 0.), 'C', mol['A'][5])),
               TesterMethod(scope='preserve_cache_del_atoms',
                            method='strip_atoms',
                            args=('H',)),
               TesterMethod(scope='preserve_cache_renumber_atoms',
                            method='renumber_atoms')
               ]

    # Test all of the tester methods
    for tester in testers:
        # The molecule cache should be empty after every run
        assert len(mol.cache) ==  0

        # Create cache objects/dicts
        obj = Mock()
        setattr(obj, tester.scope, True)
        d = {tester.scope: None}
        mol.cache['obj'] = obj
        mol.cache['d'] = d

        # Objects are in the cache
        assert 'obj' in mol.cache
        assert 'd' in mol.cache

        # Try a method that invalidates the cache in scope
        if hasattr(tester, 'args'):
            getattr(mol, tester.method)(*tester.args)
        else:
            getattr(mol, tester.method)()

        # Objects are still in the cache
        assert 'obj' in mol.cache
        assert 'd' in mol.cache

        # Reset the cached object attributes
        setattr(obj, tester.scope, False)
        del d[tester.scope]

        # The same method should now invalidate the cache
        if hasattr(tester, 'args'):
            getattr(mol, tester.method)(*tester.args)
        else:
            getattr(mol, tester.method)()

        # And the object and dict are no longer in the cache.
        assert 'obj' not in mol.cache
        assert 'd' not in mol.cache

# def test_pickle(self):
#     "Tests Pickle serialization."
#     import pickle
#
#     for name in ('2KXA', '2A5M', '5CJP'):
#         mol = Molecule(name)
#         s = pickle.dumps(mol)
#         s = pickle.dumps(mol['A'])
#         s = pickle.dumps(mol['A'][23])
#         s = pickle.dumps(mol['A'][23]['N'])
