from mollib.core.topology import topology


def test_topology_format():
    """Tests that the topology dict entries are all in the correct format"""

    # The topology keys are residue names (strings) and the values dicts.
    for key1, value1 in topology.items():
        assert isinstance(key1, str)
        assert isinstance(value1, dict)

        # The residue type should have some atoms in it
        assert len(value1) > 0

        # The value dicts shoud be the topology for each residue type. The
        # keys of those are atom names (str) and the values are sets.
        for key2, value2 in value1.items():
            assert isinstance(key2, str)
            assert isinstance(value2, set)

            # The set should have other atom names (i.e. the atom is bonded to
            # other atoms
            assert len(value2) > 0

            # All of the set items should be atom names (str)
            assert all(isinstance(i, str) for i in value2)
