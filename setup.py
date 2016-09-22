from setuptools import setup, find_packages, Extension
from Cython.Build import cythonize
import numpy as np

# Get the version number and package information
exec(open('mollib/__version__.py').read())

ext_modules = [
    Extension(name="mollib.geometry",
              sources=["mollib/core/src/geometry.pyx"],
              ) ]

setup(name=__project_name__,
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
          ],
          'distutils.commands': [
            'build_data = mollib.statistics:BuildData',
          ],
      },
      ext_modules=cythonize(ext_modules,),
      # ext_modules=cythonize('mollib/core/src/geometry.pyx'),
      include_dirs=[np.get_include()],
)
