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
