###########
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

        python setup.py build_sphinx

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

**************
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

****************
Docstring Format
****************

Docstrings follow the numpy style.

Example:

.. code-block:: python

    def calc_pa_SVD(magnetic_interactions, data):
        """Calculate the best-fit Saupe matrices for the given magnetic
        interaction arrays and RDC/RACS data.

        Parameters
        ----------
        magnetic_interactions: list of dicts
            - A list of dicts, one for each molecule to be fit.
              See :class:`mollib.pa.process_molecule.Process`.

        data: dict
            - **key**: interaction labels (str)
            - **value**: :obj:`mollib.pa.RDC` or :obj:`mollib.pa.RACS` data
              values.

        Returns
        -------
        (data_pred, Saupe_components, stats): tuple
            - data_pred: dict
                 - **key**: interaction labels (str)
                 - **value**: :obj:`mollib.pa.RDC` or :obj:`mollib.pa.RACS` data
                   values.
            - Saupe_components: dict
                 - 'S_xyz': (list of arrays) The 3x1 Saupe matrix in x/y/z repr.
                 - 'Aa': (list) The degree of alignment.
                 - 'Ar': (list) The alignment rhombicity.
                 - 'Rh': (list) The rhombicity
            - stats: dict
                 - See :func:`calc_statistics`
        """