"""
Mollib is Python package to analyze and process molecular structures and data.

Mollib is a unified command-line program and Python library for the
validation, quality analysis and manipulation of molecular structures with an
emphasis on biophysical analysis. Mollib is built on a plugin framework to
easily add new tools to manipulate and analyze structures and data, which can
then be *combined* and *cross-validated*.

Mollib includes tools for:

- The processing and protonation of molecules.
- The analysis of hydrogen bonds and hydrogen bond quality compared to the 
  highest-resolution PDB structures. Hydrogen bonds are classified based on 
  their acceptor-donor residues and Ramachandran angles
- The measurement of geometries within molecules as well as the reporting and 
  classification of Ramachandran angles.
- The statistical analysis and comparison of structures to high-resolution
  crystal structures.
- The analysis of partial alignment data with residual dipolar coupling (RDC) 
  and residual anisotropic chemical shift (RACS, a.k.a RCSA) data.
"""

from setuptools import setup, find_packages, Extension
from distutils.command.build_ext import build_ext

DOCLINES = (__doc__ or '').split("\n")

# Setup numpy include directory.
# 1. The script will not add 'numpy' to the setup_requires list, if it
#    is already installed.
# 2. Once it is installed, then the numpy include dir is added to the
#    include_dirs. This include_dirs is needed to compile the cython
#    extensions.
# 3. Finally, the install crashes if scipy is listed in the setup_requires.
#    It should be in the install_requires list, however.
setup_requires = []
try:
    import numpy as np
    include_dirs = [np.get_include(), ]
except ImportError:
    # Numpy is not available. It has to be added to the setup_requires
    setup_requires.append('numpy')
    include_dirs = []


# Setup Cython, if available
try:
    from Cython.Build import cythonize
    use_cython = True
except ImportError:
    use_cython = False


# Setup the c and c++ extensions, either with or without cython
if use_cython:
    extensions = cythonize('**/*.pyx')
else:
    import glob
    c_files = [filename for filename in glob.iglob('mollib/**/*.c')]
    extensions = []
    for c_file in c_files:
        extensions.append(Extension(c_file.strip('.c').replace('/', '.'),
                                    [c_file, ]))


# Setup a custom builder to add compilation include directories
class custom_build_ext(build_ext):

    def run(self):
        # Import numpy here, only when headers are needed
        import numpy

        # Add numpy headers to include_dirs
        self.include_dirs.append(numpy.get_include())

        # Call original build_ext command
        build_ext.run(self)


# Get the version number and package information.
# The __version__.py file is executed so that the mollib package is not loaded.
# At this point, the C/C++ extensions may not be built, and loading the mollib
# package will lead to an ImportError. This approach circumvents this problem.
__version__ = None  # This is a version string
VERSION = None  # This is a 5-item version tuple
exec(open("./mollib/__version__.py").read())


# Organize classifiers
if VERSION[3] == 'alpha':
    classifiers = ['Development Status :: 3 - Alpha', ]
elif VERSION[3] == 'beta':
    classifiers = ['Development Status :: 4 - Beta', ]
elif VERSION[3] == 'rc':
    classifiers = ['Development Status :: 4 - Beta', ]
elif VERSION[3] == 'final':
    classifiers = ['Development Status :: 5 - Production/Stable', ]
else:
    classifiers = []


classifiers += [
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
          'Operating System :: MacOS :: MacOS X',
          'Operating System :: POSIX :: Linux',
          'Programming Language :: C',
          'Programming Language :: Python',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3.5',
          'Programming Language :: Python :: 3.6',
          'Topic :: Scientific/Engineering',
          'Topic :: Scientific/Engineering :: Chemistry',
          ]


setup(name='mollib',
      version=__version__,
      license='GPLv3',
      author='Justin L Lorieau',
      description=DOCLINES[0],
      long_description="\n".join(DOCLINES[2:]),
      url='http://mollib.readthedocs.io/en/latest/',
      download_url='https://pypi.python.org/pypi/mollib/',
      include_package_data=True,
      packages=find_packages(),
      platforms='any',
      setup_requires=setup_requires,
      install_requires=['numpy', 'scipy', 'configparser'],
      tests_require=['pytest',],
      scripts=['bin/mollib', 'bin/ml'],
      entry_points={
          'console_scripts': [
            'mollib = mollib.__main__:main'
          ],
          'distutils.commands': [
            'build_data = mollib.statistics:BuildData',
          ],
      },
      ext_modules=extensions,
      include_dirs=include_dirs,
      cmdclass={'build_ext': custom_build_ext},
      zip_safe=False,
      classifiers=classifiers,
      )
