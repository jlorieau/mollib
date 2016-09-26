***********
Core Module
***********

Classes
=======

.. automodule:: mollib

.. autoclass:: Atom
    :members:

.. autoclass:: Residue
    :members:

.. autoclass:: Chain
    :members:

.. autoclass:: Molecule
    :members:

Functions
=========

.. automodule:: mollib.core

.. autofunction:: calc_vector

.. autofunction:: vector_length

.. autofunction:: within_distance

.. autofunction:: measure_distance

.. autofunction:: measure_distances

.. autofunction:: measure_angle

.. autofunction:: measure_angles

.. autofunction:: measure_dihedral

.. autofunction:: clear_cache

Plugins
=======

.. automodule:: mollib.core.plugins

    .. autoclass:: Process
        :members:

    .. autoclass:: Measure
        :members:

Settings
========

Most submodules have a ``settings.py`` file. These contain a listing of
annotated variables and values.

The SettingsManager is called by the main program (__main__) to read in
configuration files.

.. note:: Settings should be retrieve directly from the corresponding settings
          modules and not copied locally in the code. This ensures that
          settings read in from the configuration files are properly loaded
          by modules.

.. automodule:: mollib.core.settings
    :members:
