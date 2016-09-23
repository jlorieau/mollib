.. _measure-command:

Measure Command
===============
The ``measure`` command is used for measuring geometries in molecules.
All of the options and preprocessors available from the :ref:`process-command`
are also available.

    .. include:: output/cli_measure_help.html

Atom Selectors
--------------

.. _atom-selectors:

Abbreviated Selectors
~~~~~~~~~~~~~~~~~~~~~

    The measure methods find atoms using atom locators. Atom locators must
    follow one of these conventions:

        1. (residue number)-(atom name). ex: ``31-CB`` for the ``CB`` atom of
           residue number 31.
        2. (chain id).(residue number)-(atom name). ex: ``A.31-CB`` for the
           ``CB`` atom of residue number 31 in chain 'A'.

    Additionally, the chain id, residue number or both can be expressed as a
    range using the ``:`` character:

        1. (residue range)-(atom name). ex: ``31:34-CB`` for the ``CB`` atom of
           residue number 31, 32, 33 and 34.
        2. (chain range).(residue number)-(atom name). ex:``A:C.34-CB`` for the
           ``CB`` atom of residue number 34 for chains 'A', 'B', 'C' and 'D'.

    Finally, heteroatom chains have an asterisk appended to them. ex: 'C*'


    .. note:: Atom selections may encompass hundreds of atoms, which when used
              in combination, could lead to searchers over millions of
              combinations. To help improve their performance, you can either
              narrow their scope by reducing the range of chains or residue
              numbers, combine multiple _ref::`atom-filters` or use one of
              the shortcut selectors, like ``--rama`` for Ramachandran
              dihedral angles.

.. _atom-filters:

Filters
~~~~~~~

    ``--only-intra``
        Exclude atom selections that are *not* within the same residue.

    ``--exclude-intra``
        Exclude atom selections that are within the same residue.

    ``--only-intra-chain``
        Exclude atom selections that are *not* within the same chain.

    ``--exclude-intra-chain``
        Exclude atom selections that are within the same chain.

    ``--only-delta`` ``DELTA``
        Exclude atom selections that don't have at least one set of atoms
        with residues separated by ``DELTA`` number. This filter ignores the
        chain identifier and and may need to be combined
        with ``--filter-intra-chain`` or ``--exclude-intra-chain``.

    ``--only-bonded``
        Exclude atom selections that are not bonded. The bonded tests linear
        bonding relationships. For example, a dihedral with four atoms (atom1,
        atom2, atom3 and atom4) must have bonds between atom1--atom2,
        atom2--atom3 and atom3--atom4. Other bonds don't count.)

        .. note:: Bonded searches have to investigate the topology of each atom
                  selection, which can be slower than the above filters.
                  Combining the ``--only-bonded`` filter with other filters,
                  like ``--only-delta 1``, can significantly speed up searches.


Arguments
---------

    ``-d`` ``atom`` ``atom`` / ``--dist`` ``atom`` ``atom``
        Measure the distance (in Angstroms) between two atoms.

        Multiple atom pairs can used. ex: ``-d 31-N 31-CA -d 32-N 33-CA``

        Atoms must follow the standard naming conventions.
        See :ref:`atom-selectors` and :ref:`atom-filters`.

        **Examples:**

        Measure :math:`\alpha`-helical HA-H distances in chain 'A' for
        residues 23-49 of 2MUV, the homotetrametic influenza M2 channel. Include
        statistics on the measured distances.

        .. include:: output/cli_measure_i_2MUV_d_23:49-HA_23:49-H_only-delta_3_stats.html

        Measure CA-CA distances between residue 20-21 for chains 'A', 'B', 'C'
        and 'D' of 2MUV--excluding same residue distances and same chain
        distances.

        .. include:: output/cli_measure_i_2MUV_d_A:D.20:21-CA_A:D.20:21-CA_exclude-intra_exclude-intra-chain.html

        Compare the distance between the HA of residue 5 and the H of residue
        21 for two different structures, 2KXA and 2LWA. The 2KXA structure
        represents the wildtype hemagglutinin fusion peptide (HAfp) in the
        *closed* helical-hairpin structure, placing these two atoms in close
        promixity. The 2LWA structure represents the conformational ensemble
        of theHAfp-G8A mutant with a closed structure (chain 'A'), and
        semi-closed structure (chain 'B') and an open structure (chain 'C').

        .. include:: output/cli_measure_i_2KXA_2LWA_d_A:C.5-HA_A:C.21-H_only-intra-chain.html

    ``-a`` / ``--angle``
        Measure the angle (in degrees) between three atoms.

        Multiple atom triplets can be used. ex: ``-a 31-N 31-CA 31-CB
        -a 32-N 32-CA 32-CB``

        Atoms must follow the standard naming conventions.
        See :ref:`atom-selectors` and :ref:`atom-filters`.

        **Examples:**

        Measure the angle of the bonded 'C-1'--'N'--'H' atoms for residues
        20-30 from the ubiquitin structure 2MJB.

        .. include:: output/cli_measure_i_2MJB_a_20:30-C_20:30-N_20:30-H_only-bonded.html

    ``-dih`` / ``--dihedral``
        Measure the dihedral angle (in degrees) between four atoms.

        Multiple atom quartets can be used. ex: ``-dih 30-C 31-N 31-CA 31-C
        -dih 31-N 31-CA 31-C 32-N``

        Atoms must follow the standard naming conventions.
        See :ref:`atom-selectors` and :ref:`atom-filters`.


        .. note:: If simple Ramachandran and side-chain
                  dihedrals are needed, checkout ``--rama``, ``--chi-1``.


        **Examples**

        .. include:: output/cli_measure_i_2KXA_dih_2:6-C_2:6-N_2:6-CA_2:6-C_only-bonded_stats.html

    ``--rama``
        Measure Ramachandran angles (in degrees) for a protein. Filters and
        options are ignored. Heteroatom chains are skipped.

        The ``--rama`` command classifies Ramachandran angles based on
        backbone-backbone amide hydrogen bonds. A residue is classified based
        on whether its amide or carbonyl is participating in a hydrogen bond.
        Residues without a classification are either randomly coil, or they
        correspond to secondary structure units at the surface of the protein,
        without an intramolecular hydrogen bond.

        The *isolated* classification is given for residues that have backbone
        hydrogen bonds, but these cannot be classified into conventional
        secondary structure units. See the :ref:`hbonds_command` for further
        details.

        **Examples**

        Measure the Ramachandran :math:`\phi` and :math:`\psi` angles for the
        hemagglutinin fusion peptide structure 2KXA.

        .. include:: output/cli_measure_i_2KXA_rama.html

Options
~~~~~~~

    ``--stats``
        Report the average and standard deviation of all measured values.