from setuptools import setup, find_packages

import mollib

setup(name='mollib',
      license='GPL V2',
      author='Justin L Lorieau',
      description='A Python package to read molecular structures',
      packages=['mollib'],
      platforms='any',
      test_suite='nose.collector',
      entry_points={
          'console_scripts': [
              'mollib = mollib.__main__:main'
          ]
      },
)
