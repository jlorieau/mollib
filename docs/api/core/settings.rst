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