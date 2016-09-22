from setuptools import setup, find_packages, Extension
from Cython.Build import cythonize
import numpy as np

__name__ = 'mollib'
__author__ = 'Justin L Lorieau'
__versioninfo__ = (1, 3, 'a1')
__version__ = '.'.join(map(str, __versioninfo__))

ext_modules = [
    Extension(name="mollib.geometry",
              sources=["mollib/core/src/geometry.pyx"],
              ) ]

setup(name=__name__,
      version=__version__,
      license='GPL V2',
      author= __author__,
      description='A Python package to read molecular structures',
      install_requires=['numpy'],
      include_package_data = True,
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
