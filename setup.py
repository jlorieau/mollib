from setuptools import setup, find_packages, Extension
from Cython.Build import cythonize
import numpy as np

# import mollib

# find . -name "*.so" | xargs rm
# python setup.py build_ext --inplace  # python2
# python3 setup.py build_ext --inplace  # python3

ext_modules = [
    Extension(name="mollib.geometry",
              sources=["mollib/core/src/geometry.pyx"],

              ) ]

setup(name='mollib',
      # version=mollib.__version__,
      license='GPL V2',
      author='Justin L Lorieau',
      description='A Python package to read molecular structures',
      packages=find_packages(),
      platforms='any',
      test_suite='nose.collector',
      scripts=['bin/mollib'],
      entry_points={
          'console_scripts': [
              'mollib = mollib.__main__:main'
          ]
      },
      ext_modules=cythonize(ext_modules,),
      # ext_modules=cythonize('mollib/core/src/geometry.pyx'),
      include_dirs=[np.get_include()],
)
