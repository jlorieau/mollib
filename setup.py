from setuptools import setup, find_packages, Extension
import numpy as np


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


# Get the version number and package information.
# The __version__.py file is executed so that the mollib package is not loaded.
# At this point, the C/C++ extensions may not be built, and loading the mollib
# package will lead to an ImportError. This approach circumvents this problem.
__version__ = None  # This is a version string
VERSION = None  # This is a 5-item version tuple
execfile('mollib/__version__.py')  # The following loads __version__/VERSION


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
          'Operating System :: OS Independent',
          'Programming Language :: C',
          'Programming Language :: Python',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3.4',
          'Programming Language :: Python :: 3.5',
          'Programming Language :: Python :: 3.6',
          'Topic :: Scientific/Engineering',
          'Topic :: Scientific/Engineering :: Chemistry',
          ]


setup(name='mollib',
      version=__version__,
      license='GPLv3',
      author='Justin L Lorieau',
      description=('A Python package to analyze and process molecular '
                   'structures and data'),
      include_package_data=True,
      packages=find_packages(),
      platforms='any',
      install_requires=['numpy', 'scipy'],
      tests_requires=['nose>=1.0'],
      test_suite='nose.collector',
      scripts=['bin/mollib'],
      entry_points={
          'console_scripts': [
            'mollib = mollib.__main__:main'
          ],
          'distutils.commands': [
            'build_data = mollib.statistics:BuildData',
          ],
      },
      ext_modules=extensions,
      include_dirs=[np.get_include()],
      classifiers=classifiers,
      )
