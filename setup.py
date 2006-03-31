#!/usr/bin/env python


import distutils
from distutils.core import setup, Extension

from  tdl.version import name,version,author,email,desc
setup(
    name = name,
    version = version,
    author =  author,
    author_email = email,
    description = desc,
    package_dir = {'tdl': 'tdl'},
    packages = ['tdl'],
    license = 'Python',
)
