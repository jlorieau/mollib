Development
###########

The following are instructions for locally developing a branch of mollib.

1. **Checkout**. First, checkout the branch using git.

2. **Build Cython Extensions**. The Cython and C extensions can be compiled in
   place.

    ..  code-block:: shell-session

        python setup.py build_ext --inplace

3. **Build Datasets** Build and compute the datasets. This may take a few hours.

    .. code-block:: shell-session

        python setup.py build_data

4. **Build Documentation** Build the documentation in html format under
   ``docs/``

    .. code-block:: shell-session

        make docs

    .. note:: You will need to have ``sphinx``, ``sphinxcontrib-napoleon``
              installed.

5. **Test the Build** Run the package's tests.

    .. code-block:: shell-session

        python setup.py nosetests

    .. note:: You will need ``nosetests`` to run the tests.

6. **Install in Developer Mode**. Install the package in developer mode. This
adds the package's source path to the python path. Edits to the source path
are reflected in the global script.

    .. code-block:: shell-session

        python setup.py develop

    If you'd like to uninstall the develop mode, use the following command.

    .. code-block:: shell-session

        python setup.py develop --uninstall

Building Platform Packages
**************************

Platform specific distributions are needed because mollib includes C extensions
written in Cython. These extension must be compiled using a compiler for each
type of operating system.

The build and install distribution packages, the following packages are needed:

    - `Cython 0.25+ <http://cython.org>`_

Linux Distributions (Python Wheels)
===================================

    Compile and build a Python wheel package.

    .. code-block:: shell-session

        python setup.py bdist_wheel

    The Python wheel can be installed using
    `pip <https://pypi.python.org/pypi/pip>`_.

    .. code-block:: shell-session

        sudo pip install <package_file.whl>

Mac OS X Distributions (mpkg)
=============================

    Compile and build a Mac OS X package (``.mpkg``) file.

    .. code-block:: shell-session

        python setup.py bdist_mpkg

    The ``.mpkg`` file can be installed by dragging it into the Applications
    folder.


Adding Plugins
**************

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
****************

Docstrings follow the numpy style.

Example:

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