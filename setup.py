'''
Setup file for package snapperdb

:author: alex
:created: 18 Dec 2014
'''

from distutils.core import setup

setup(name='snapperdb',
      version='1.0',
      description='Description goes here',
      author='Author information',
      author_email='email',
      url='URL of the repo?',
      packages=['snapperdb', 'snapperdb.gbru_vcf', 'snapperdb.snpdb'],
      requires=["psycopg2"]  # This needs to be extended to include others.
     )