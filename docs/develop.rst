Development
###########

The following are instructions for locally developing a branch of mollib.

1. **Checkout**. First, checkout the branch using git.

2. **Build Cython Extensions**. The Cython and C extensions can be compiled in
   place.

    ..  code-block:: shell-session

        $ make inplace

3. **Build Datasets** Build and compute the datasets. This may take a few hours.

    .. code-block:: shell-session

        $ make build-data

4. **Build Documentation** Build the documentation in html format under
   ``docs/``

    .. code-block:: shell-session

        $ make docs

    .. note:: You will need to have ``sphinx``, ``sphinxcontrib-napoleon``,
              ``latex`` and ``latexmk`` installed.

5. **Test the Build** Run the package's tests.

    .. code-block:: shell-session

        $ make test

    .. note:: You will need ``nosetests`` to run the tests.

    Alternatively, mollib can be tested against multiple platforms using
    `tox <https://tox.readthedocs.io/en/latest/>`_.

    .. code-block:: shell-session

        $ make test-all

6. **Install in Developer Mode**. Install the package in developer mode. This
   adds the package's source path to the python path. Edits to the source path
   are reflected in the global script.

    .. code-block:: shell-session

        $ make develop

    If you'd like to uninstall the develop mode, use the following command.

    .. code-block:: shell-session

        $ python setup.py develop --uninstall

Testing
=======

Mollib includes 4 different kinds of tests. These are all executed by the
``make test-all`` command.

1. **Pytests and unitests**. These are tests stored in the ``tests`` directory.
   They are run automatically by executing one of the following:

        a. ``pytest``
        b. ``make test``
        c. ``make test-all``

2. **Docstring tests**. These are tests within the docstrings of functions,
   classes and methods in the mollib source (``mollib`` subdirectory). These
   tests are run automatically by executing one of the following:

        a. ``pytest``
        b. ``make test``
        c. ``make test-all``

3. **Tox**. All of the pytests, unittests and doctests are tested in each Python
   environment using tox. These tests are run automatically by executing one of
   the following:

        a. ``tox``
        b. ``make test-all``

4. **CLI tests**. These tests detect changes in the output text from a specific
   set of command arguments and mollib. These tests are located in
   ``tests/cli`` directory. An error will be flagged if the output of a mollib
   command has changed from the contents of the ``.txt`` file. The output of
   the commands are also used in the documentation within the ``.rst`` files.
   These tests are run by executing:

        a. ``make test-cli``
        b. Additionally, if the output of a command changes, the reference
           commands can be reset by entering the ``tests/cli`` directory and
           typing ``make clean&&make build``. The changed ``.txt`` and ``.rst``
           files should be committed to the repository.
        c. New commands can be created by typing the command in a ``.sh`` file
           and making this file executable.

Makefile Options
================

The ``make`` command contains a number of commands to setup and develop
mollib. The make commands are listed by typing ``make help``.

Including Datasets
==================

Datasets are included in the ``mollib/data`` directory. Data files should be
included in the ``MANIFEST.in`` file.

Adding Plugins
==============

Adding plugin modules may require the registration of the plugin, the
registration of the plugin's settings or both.

To register a plugin, add the following to the root ``__init__.py`` file for
the plugin:

    .. code-block:: python

        from .plugin import PluginClass
        plugin = PluginClass()

To register a plugin's settings, add the following to the root ``__init__.py``
file for the plugin:

    .. code-block:: python

        from . import settings

        from mollib.core import register_settings
        register_settings(settings)

Branches
========

Mollib uses git-flow to organize branches. These are the main branches:

1. **master**. The master branch is the default branch, and it contains
   the *production* code. It should pass all tests including travis-ci tests.

2. **development**. The development branch is used to accept new features,
   hotfixes and releases. It should pass all tests including travis-ci tests.

Docstring Format
================

Docstrings follow the numpy style. There are a few additional guidelines:

    1. ``dict`` parameters and return values should list the expected
       keys/values

      a. ``dict`` parameters should list the key and value types using 'key' and
         'value' in bold. If known, the object type should be listed after the
         description.

        .. code:: raw

            - **key**: interaction label, str

      b. ``dict`` return values should either list the key/value pairs, or list
         specific keys and values.

        .. code:: raw

            - 'Q (%)': The fit Q-factor in percentage, float

    2. *Sublists* should have a new line before the sublisting.

        .. code:: raw

           - 'Overall': Overall Statistics, :obj:`collections.OrderedDict`

              - 'Q (%)': The fit Q-factor in percentage, float
              - 'RMS': The root-mean square of the fit (Hz/ppb), float
              - 'count': The number of interactions fit, int

    3. *Lists* that follow a paragraph listing in a parameter should not be
       indented with respect to the paragraph.

        .. code:: raw

            angles: dict
                A dict of the angles between atoms that define the hydrogen
                bond.

                - **key**: tuple of three :obj:`Atom` objects
                - **value**: the angle (in deg) between the :obj:`Atom` objects

    4. *args* and *kwargs* args are listed separately and as optional
       parameters.

        .. code:: raw

            Parameters
            ----------
            args: tuple, optional
                If specified a default argument, then this will be returned if
                the key isn't found. Otherwise a ValueError exception is raised.
            kwargs: dict, optional
                If specified a default argument, then this will be returned if
                the key isn't found. Otherwise a ValueError exception is raised.

Example 1
*********

.. code-block:: python

   def calc_summary(magnetic_interactions, Saupe_components, data, predicted):
        """Calculate the statistics between predicted and calculated RDCs and
        RACSs.

        Parameters
        ----------
        magnetic_interactions: list of dicts
            - A list of dicts, one for each molecule to be fit.
              See :class:`mollib.pa.process_molecule.Process`
        Saupe_components: dict
            See the output of :func:`mollib.pa.svd.calc_pa_SVD`
        data: dict
            - **key**: interaction labels, str
            - **value**: :obj:`mollib.pa.RDC` or :obj:`mollib.pa.RACS` data
              values.
        predicted: dict
            - **key**: interaction labels, str
            - **value**: :obj:`mollib.pa.RDC` or :obj:`mollib.pa.RACS` data
              values.

        Returns
        -------
        summary: :obj:`collections.OrderedDict`

            - 'Overall': Overall Statistics, :obj:`collections.OrderedDict`

              - 'Q (%)': The fit Q-factor in percentage, float
              - 'RMS': The root-mean square of the fit (Hz/ppb), float
              - 'count': The number of interactions fit, int

            - 'Alignment': Details on the alignment tensor,
              :obj:`collections.OrderedDict`

              - 'Aa': The alignment tensor anisotropy, float
              - 'Ar': The alignment tensor rhobicity, float

            - 'Saupe': Details on the Saupe matrix, :obj:`collections.OrderedDict`

              - 'Szz': The zz-component of the Saupe matrix, float
              - 'Sxx': The xx-component of the Saupe matrix, float
              - 'Syy': The yy-component of the Saupe matrix, float

            - 'Angles': Alignment tensor orientation in Rose convention,
              :obj:`collections.OrderedDict`

              - "Z (deg)": The alignment alpha angle (deg), float
              - "Y' (deg)": The alignment beta angle (deg), float
              - "Z'' (deg)": The alignment gamma angle (deg), float
        """

Example 2
*********

.. code-block:: python

    def fill_gaps(molecule, classifications, classification_type, dihedral_test,
                  extend_terminii=False, label_N_term=0, label_C_term=0,
                  gap_tolerance=1, overwrite_assignments=False):
        """Fill gaps in the classifications dict assignments.

        Gaps occur in the secondary structure assignment from hydrogen bonds,
        for example, with beta-strands on the edges of beta sheets. This
        function finds stretches of secondary structure assignments, it checks
        the dihedral angles and fills in gaps in the stretches. For a sheet:
        'E E E E' becomes 'EEEEEEE'.

        Parameters
        ----------
        molecule: :obj:`mollib.Molecule`
            The molecule object to classify the secondary structure elements.
        classifications: dict
            A dict with the classifications.

              - **key**: (chain.id, residue.number). ex: ('A', 31)
              - **value**: (major_classification, minor_classification).
                ex: ('alpha-helix', 'N-term')
        classification_type: str
            The name of the classification type to check. ex: 'alpha-helix'
        dihedral_test: function or None
            - A test function that takes a :obj:`mollib.Residue` and returns
              True if the residue's dihedral angles are within range for the
              'classification_type'.
            - If None is specified, then the dihedral angles of residues will
              not be tested.
        extend_terminii: bool or int, optional
            If True, the previous and subsequence residues of each contiguous
            stretch of residue classification will be checked to see if they fall
            within the dihedral angle range as well.
        label_N_term: int (optional)
            Label the first 'N' number of residues in the contiguous block as
            'N-term'
        label_C_term: int, optional
            Label the last 'N' number of residues in the contiguous block as
            'C-term'
        gap_tolerance: int, optional
            The assignment of contiguous stretches of a secondary structure
            assignment will tolerate this number of 'gaps' in the residue
            numbers.
            For a gap_toleranace of 1 and a checked sheet assignment, the
            following group 'E E E E' will be treated as a single contiguous
            block of sheetassignments.
        overwrite_assignments: bool, optional
            If True, classification assignments will be overwritten, if an
            assignments has already been made for a given residue.

        Returns
        -------
        None
        """