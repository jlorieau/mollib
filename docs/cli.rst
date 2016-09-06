=====================
Commandline Interface
=====================
The command line interface includes all of the mollib functions for processing
molecules.

    .. literalinclude:: cli/cli_help.txt
        :language: shell-session

    ``-h`` / ``--help``
        Get basic help information on usage of the command line interface, or one
        one of the mollib commands.

    ``-d`` / ``--debug``
        Display debug messages to the terminal. This is usually only useful to
        developers or in helping to debug issues.

    ``-v`` / ``--verbose``
        Display informational messages to the terminal. By default, only warning
        and error messages are displayed.

    ``--version``
        Display the version number of the installed mollib.

********
Commands
********

Commands include various types of operations that can be conducted by mollib.
The command listing may contain additional commands, depending on whether
mollib plugins are installed.

.. _process-command:

Process Command
===============
The ``process`` command is the main command for processing, reading and writing
files. All of the options and preprocessors available in ``process`` are
available to other commands.

    .. literalinclude:: cli/cli_process_help.txt
        :language: shell-session

Arguments
---------

    ``-i`` ``id/filename`` / ``--in`` ``id/filename``
        **(required)** The listing of one or more structural identifiers
        (ex: PDB file identifiers) or filenames.

        If the structure could not be found locally, a copy will be
        downloaded and cached for further analysis.

        Multiple input identifiers and filename can be used simultaneously to
        process multiple files.

    ``-o`` ``filename`` / ``--out`` ``filename``
        The output filename(s).

        Structure files that are written have passed through the mollib parser and
        will likely not be identical to the original files. Changes may have been
        done to the atoms, header comments or other aspects, depending on which
        preprocessors were used.

        Multiple output filenames can be used, and these will be matched
        to the corresponding entries in the input filenames or identifiers.

    ``-c`` ``filename`` / ``--config`` ``filename``
        The configuration file.

Preprocessors Arguments
-----------------------

    ``--hydrogenate``
        Strip all hydrogen atoms and re-add hydrogens based on ideal geometry.
        More details on the methods can be found in the API documentation
        :doc:`plugins/hydrogens`.


        .. note:: Adding two hydrogens to an sp2 heavy atom will label the
                  E-hydrogen 'H1' and the Z-hydrogen 'H2'. This situation
                  happens with the HD2 hydrogens of asparagine residues, for
                  example.

        .. note:: Adding two hydrogens to an sp3 heavy atom should will label
                  the pro-R hydrogen 'H2' and the pro-S hydrogen 'H3'. Some
                  exceptions in proteins exist. By default, the JBNMR 12, 1-23
                  (1998) convention is followed.

.. _measure-command:

Measure Command
===============
The ``measure`` command is used for measuring geometries in molecules.
All of the options and preprocessors available from the :ref:`process-command`
are also available.

    .. literalinclude:: cli/cli_measure_help.txt
        :language: shell-session

Atom Selectors
--------------

.. _atom-selectors:

Abbreviated Selectors
~~~~~~~~~~~~~~~~~~~~~

    The measure methods find atoms using atom locators. Atom locators must
    follow one of these conventions:

        1. (residue number)-(atom name). ex: ``31-CB`` for the ``CB`` atom of
           residue number 31.
        2. (chain id)-(residue number)-(atom name). ex: ``A.31-CB`` for the
           ``CB`` atom of residue number 31 in chain 'A'.

    Additionally, the chain id, residue number or both can be expressed as a
    range using the ``:`` character:

        1. (residue range)-(atom name). ex: ``31:34-CB`` for the ``CB`` atom of
           residue number 31, 32, 33 and 34.
        2. (chain range)-(residue number)-(atom name). ex:``A:C.34-CB`` for the
           ``CB`` atom of residue number 34 for chains 'A', 'B', 'C' and 'D'.

    Finally, heteroatom chains have an asterisk appended to them. ex: 'C*'

.. _atom-filters:

Filters
~~~~~~~

    ``--filter-intra``
        Exclude atom selections that are *not* within the same residue number.
        This filter ignores the chain identifier and may need to be combined
        with ``--filter-intra-chain`` or ``--exclude-intra-chain``.

    ``--exclude-intra``
        Exclude atom selections that are within the same residue number.
        This filter ignores the chain identifier and and may need to be combined
        with ``--filter-intra-chain`` or ``--exclude-intra-chain``.

    ``--filter-intra-chain``
        Exclude atom selections that are *not* within the same chain.

    ``--filter-delta`` ``DELTA``
        Exclude atom selections that don't have at least one set of atoms
        with residues separated by ``DELTA`` number. This filter ignores the
        chain identifier and and may need to be combined
        with ``--filter-intra-chain`` or ``--exclude-intra-chain``.

    ``--filter-bonded``
        Exclude atom selections that are not bonded. The bonded tests linear
        bonding relationships. For example, a dihedral with four atoms (atom1,
        atom2, atom3 and atom4) must have bonds between atom1--atom2,
        atom2--atom3 and atom3--atom4. Other bonds don't count.)

Arguments
---------

    ``--exclude-intra``
        Exclude measurements within the same residue

    ``--delta``
        Report measurements only for residues separated by DELTA number of
        residues.

    ``-d`` ``atom`` ``atom`` / ``--dist`` ``atom`` ``atom``
        Measure the distance (in Angstroms) between two atoms.

        Multiple atom pairs can used. ex: ``-d 31-N 31-CA -d 32-N 33-CA``

        Atoms must follow the standard naming conventions.
        See :ref:`atom-selectors` and :ref:`atom-filters`.

        **Examples:**

        Measure :math:`\alpha`-helical HA-H distances in chain 'A' for residues 23-49.

        .. literalinclude:: cli/cli_measure_i_2MUV_d_23:49-HA_23:49-H_only-delta_3.txt
            :language: shell-session

        Measure CA-CA distances between residue 20-21 for chains 'A', 'B', 'C'
        and 'D'--excluding same residue distances and same chain distances

        .. literalinclude:: cli/cli_measure_i_2MUV_d_A:D.20:21-CA_A:D.20:21-CA_exclude-intra_exclude-intra-chain.txt
            :language: shell-session

        Measure

        .. literalinclude:: cli/cli_measure_i_2KXA_2LWA_d_A:C.5-HA_A:C.21-H_only-intra-chain.txt
            :language: shell-session

    ``-a`` / ``--angle``
        Measure the angle (in degrees) between three atoms.

        Multiple atom triplets can be used. ex: ``-a 31-N 31-CA 31-CB
        -a 32-N 32-CA 32-CB``

        Atoms must follow the standard naming conventions.
        See :ref:`atom-selectors` and :ref:`atom-filters`.

    ``-dih`` / ``--dihedral``
        Measure the dihedral angle (in degrees) between four atoms.

        Multiple atom quartets can be used. ex: ``-dih 30-C 31-N 31-CA 31-C
        -dih 31-N 31-CA 31-C 32-N``

        Atoms must follow the standard naming conventions.
        See :ref:`atom-selectors`.

    ``-r`` / ``--ramachandran``
        Display a (Markdown) table of the structure's ramachandran angles
        (in degrees).