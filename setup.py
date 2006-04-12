#!/usr/bin/env python


import distutils
from distutils.core import setup, Extension

from  lib.version import name,version,author,email,desc
setup(
    name = name,
    version = version,
    author =  author,
    author_email = email,
    description = desc,
    package_dir = {'tdl': 'lib'},
    packages = ['tdl','tdl.TkPlotter'],
    license = 'Python',
    data_files  = [('bin',['tdl'])]
)
