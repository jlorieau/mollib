.. _process-command:

``process`` command
===================
The ``process`` command is the main command for processing, reading and writing
files. All of the options and preprocessors available in ``process`` are
available to other commands.

.. include:: cmds/ml_process_help.rst

Arguments
---------

``-i`` / ``--in`` ``id/filename``
    **(required)** One or more structural identifiers (ex: PDB file
    identifiers) or filenames to load.

    If the structure could not be found locally, a copy will be
    downloaded and cached for further analysis.

    Multiple input identifiers and filenamed can be used to process multiple
    files simultaneously.

``-o`` / ``--out`` ``filename``
    The output filename(s).

    Structure files that are written have passed through the mollib parser and
    will likely not be identical to the original files. Changes may have been
    done to the atoms, header comments or other aspects, depending on which
    preprocessors were used.

    Multiple output filenames can be used, and these will be matched
    to the corresponding entries in the input filenames or identifiers.

``-c`` / ``--config`` ``filename``
    The configuration file. See :ref:`Configuration
    files <configuration-files>` for more details.

``-s`` / ``--save``
    Save fetched files from the internet to the local directory.

``-m`` / ``--models``
    The specific model numbers to load. By default, only the first model is
    loaded.

``-l``
    List details on the molecule, including the number of chains, residues
    and atoms. This option is helpful in determining selections in
    :ref:`atom-selectors`.

Preprocessors Arguments
^^^^^^^^^^^^^^^^^^^^^^^

.. _hydrogenate:

``--hydrogenate``
    Strip all hydrogen atoms and re-add hydrogens based on ideal geometry.
    More details on the current implementation can be found in the API
    documentation :doc:`../api/hydrogens/hydrogens`.


    .. note:: Adding two hydrogens to an sp2 heavy atom will label the
              E-hydrogen 'H1' and the Z-hydrogen 'H2'. This situation
              happens with the HD2 hydrogens of asparagine residues, for
              example.

    .. note:: Adding two hydrogens to an sp3 heavy atom will label
              the pro-R hydrogen 'H2' and the pro-S hydrogen 'H3'. Some
              exceptions in proteins exist. By default, the JBNMR 12, 1-23
              (1998) convention is followed.

Examples
--------

The following example loads a crystal structure of ubiquitin from the PDB
(``-i 1UBQ``), adds hydrogens to the molecule (``--hydrogenate``) and saves the
output to a new file (``-o 1UBQ_H.pdb``).

.. include:: cmds/ml_process_1ubq_1.rst
