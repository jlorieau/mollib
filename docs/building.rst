Building and Deployment
#######################

Building
========

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

    1. Upload the source distribution

        .. code-block:: shell-session

            $ python setup.py sdist upload

    2. Upload a binary distribution

        .. code-block:: shell-session

            $ python setup.py bdist_wheel upload

    2. Compile and build a Mac OS X package (``.mpkg``) file.

        .. code-block:: shell-session

            $ python setup.py bdist_mpkg

      The ``.mpkg`` file can be installed by dragging it into the Applications
      folder.

Deployment
==========

The linux and osx packages are built using the mollib-wheels repository. Follow
these steps to create a new release.

1. **mollib repository**. Tag a new release on the mollib *master* branch and
   bump the version in the ``mollib/__version__.py`` file.

2. **mollib-wheelhouse**. Update the ``BUILD_COMMIT`` in the ``.travis.yml`` to
   point to the new mollib tag. If needed, switch the twine repository to the test
   servers.

3. Push the sdist to PyPi.
