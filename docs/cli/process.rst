.. _process-command:

Process Command
===============
The ``process`` command is the main command for processing, reading and writing
files. All of the options and preprocessors available in ``process`` are
available to other commands.

    .. literalinclude:: output/cli_process_help.txt
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
        :doc:`../plugins/hydrogens`.


        .. note:: Adding two hydrogens to an sp2 heavy atom will label the
                  E-hydrogen 'H1' and the Z-hydrogen 'H2'. This situation
                  happens with the HD2 hydrogens of asparagine residues, for
                  example.

        .. note:: Adding two hydrogens to an sp3 heavy atom should will label
                  the pro-R hydrogen 'H2' and the pro-S hydrogen 'H3'. Some
                  exceptions in proteins exist. By default, the JBNMR 12, 1-23
                  (1998) convention is followed.