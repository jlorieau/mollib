from setuptools import setup, find_packages

import mollib

setup(name='mollib',
      version=mollib.__version__,
      license='GPL V2',
      author='Justin L Lorieau',
      description='A Python package to read molecular structures',
      packages=find_packages(),
      platforms='any',
      test_suite='nose.collector',
      entry_points={
          'console_scripts': [
              'mollib = mollib.__main__:main'
          ]
      },
)
