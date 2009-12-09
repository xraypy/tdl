#!/usr/bin/env python


import distutils
from distutils.core import setup, Extension


setup(
    name = 'larch',
    version = '0.8.0',
    author = 'Matthew Newville',
    author_email = 'newville@cars.uchicago.edu',
    license = 'Python',
    description = 'A data processing macro language for python',
    package_dir = {'larch': 'lib',
                   'larch.modules':'modules',
                   'larch.modules':'modules',
                   },
    packages = ['larch','larch.modules'],
    package_data = {'larch.modules':['startup.lar']},
    data_files  = [('bin',['larch'])],)

