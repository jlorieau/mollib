###########
Development
###########

The following are instructions for locally developing a branch of mollib.

First, checkout the branch using git.

The Cython and C extensions can be compiled in place.

..  code-block:: shell-session

    python setup.py build_ext --inplace

Install the package in developed mode. This adds the package's source path
to the python path. Edits to the source path are reflected in the global
script.

.. code-block:: shell-session

    python setup.py develop

If you'd like to uninstall the develop mode, use the following command.

.. code-block:: shell-session

    python setup.py develop --uninstall

To test the package, run the nosetests

.. code-block:: shell-session

    python setup.py nosetests

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