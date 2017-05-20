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

    .. note:: You will need to have ``sphinx``, ``sphinxcontrib-napoleon``
              installed.

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

Makefile Options
================

The ``make`` command contains a number of commands to setup and develop
mollib.

    .. include:: cli/output/make_help.rst

Building Platform Packages
==========================

Platform specific distributions are needed because mollib includes C extensions
written in Cython. These extension must be compiled using a compiler for each
type of operating system.

The build and install distribution packages, the following packages are needed:

    - `Cython 0.25+ <http://cython.org>`_

Linux Distributions (Python Wheels)
***********************************

    Compile and build a Python wheel package.

    .. code-block:: shell-session

        $ python setup.py bdist_wheel

    The Python wheel can be installed using
    `pip <https://pypi.python.org/pypi/pip>`_.

    .. code-block:: shell-session

        $ sudo pip install <package_file.whl>

Mac OS X Distributions (mpkg)
*****************************

    Compile and build a Mac OS X package (``.mpkg``) file.

    .. code-block:: shell-session

        $ python setup.py bdist_mpkg

    The ``.mpkg`` file can be installed by dragging it into the Applications
    folder.


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

Docstring Format
================

Docstrings follow the numpy style.

Example 1
*********

.. code-block:: python

    def calc_statistics(magnetic_interactions, Saupe_components, data, predicted):
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
            - **key**: interaction labels (str)
            - **value**: :obj:`mollib.pa.RDC` or :obj:`mollib.pa.RACS` data
              values.
        predicted: dict
            - **key**: interaction labels (str)
            - **value**: :obj:`mollib.pa.RDC` or :obj:`mollib.pa.RACS` data
            values.

        Returns
        -------
        stats: :obj:`collections.OrderedDict`
            - 'Q': (float) the Q-factor of the fit
            - 'R': (float) the R-factor of the fit
            - 'RMS': (Hz/ppb) the root-mean square of the fit
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
        extend_terminii: bool or int (optional)
            If True, the previous and subsequence residues of each contiguous
            stretch of residue classification will be checked to see if they fall
            within the dihedral angle range as well.
        label_N_term: int (optional)
            Label the first 'N' number of residues in the contiguous block as
            'N-term'
        label_C_term: int (optional)
            Label the last 'N' number of residues in the contiguous block as
            'C-term'
        gap_tolerance: int (optional)
            The assignment of contiguous stretches of a secondary structure
            assignment will tolerate this number of 'gaps' in the residue
            numbers.
            For a gap_toleranace of 1 and a checked sheet assignment, the
            following group 'E E E E' will be treated as a single contiguous
            block of sheetassignments.
        overwrite_assignments: bool (optional)
            If True, classification assignments will be overwritten, if an
            assignments has already been made for a given residue.

        Returns
        -------
        None
        """