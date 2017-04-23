from setuptools import setup, find_packages, Extension
import numpy as np

# Get the version number and package information
version = __import__('mollib').__version__

# Setup Cython, if available
try:
    from Cython.Build import cythonize
    use_cython = True
except ImportError:
    use_cython = False

ext = '.pyx' if use_cython else '.c'

# extensions = [
#     Extension(["mollib/core/geometry" + ext],
#               ) ]

if use_cython:
    extensions = cythonize('**/*.pyx')


setup(name='mollib',
      version=version,
      license='GPL V2',
      author='Justin L Lorieau',
      description='A Python package to read molecular structures',
      install_requires=['numpy', 'scipy'],
      include_package_data = True,
      packages=find_packages(),
      platforms='any',
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
      # ext_modules=cythonize("**/*.pyx"),
      include_dirs=[np.get_include()],
)
