Overview
========
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

Commands
========

Commands include various types of operations that can be conducted by mollib.
The command listing may contain additional commands, depending on whether
mollib plugins are installed.

.. _process-command:

Process Command
---------------
The ``process`` command is the main command for processing, reading and writing
files. All of options and preprocessors are available in ``process`` and other
commands.

    .. literalinclude:: cli/cli_process_help.txt
        :language: shell-session

Basic Arguments
~~~~~~~~~~~~~~~

    ``-i`` / ``--in``
        The listing of one or more structural identifiers (ex: PDB file identiers)
        or filenames.

        If the structure could not be found locally, a copy will be
        downloaded and cached for further analysis.

        Multiple input identifiers and filename can be used simultaneously to
        process multiple files.

    ``-o`` / ``--out``
        The output filename(s).

        Structure files that are written have passed through the mollib parser and
        will likely not be identical to the original files. Changes may have been
        done to the atoms, header comments or other aspects, depending on which
        preprocessors were used.

        Note that multiple output filenames can be used, and these will be matched
        to the corresponding entries in the input filenames or identifiers.

    ``-c`` / ``--config``
        The configuration file.

Preprocessors Arguments
~~~~~~~~~~~~~~~~~~~~~~~

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
---------------
The ``measure`` command is used for measuring geometries in molecules.
All of options and preprocessors are available in the :ref:`process-command`
are also available.

    .. literalinclude:: cli/cli_measure_help.txt
        :language: shell-session

Basic Arguments
~~~~~~~~~~~~~~~

    ``-d`` / ``--dist``
        Measure the distance (in Angstroms) between two atoms.

        Multiple atom pairs can used. ex: ``-d 31-N 31-CA -d 32-N 33-CA``

        .. _atom-conventions:

        Atoms must follow one of the following *atom conventions*:
            1. (residue number)-(atom name). ex: ``31-CB``
            2. (chain number)-(residue number)-(atom name). ex: ``A-31-CB``

    ``-a`` / ``--angle``
        Measure the angle (in degrees) between three atoms.

        Multiple atom triplets can be used. ex: ``-a 31-N 31-CA 31-CB
        32-N 32-CA 32-CB``

        Atoms must follow the standard naming conventions.
        See :ref:`atom conventions <atom-conventions>`.

    ``-dih`` / ``--dihedral``
        Measure the dihedral angle (in degrees) between four atoms.

        Multiple atom quartets can be used. ex: ``-dih 30-C 31-N 31-CA 31-C
        -dih 31-N 31-CA 31-C 32-N``

        Atoms must follow the standard naming conventions.
        See :ref:`atom conventions <atom-conventions>`.

    ``-r`` / ``--ramachandran``
        Display a (Markdown) table of the structure's ramachandran angles
        (in degrees).