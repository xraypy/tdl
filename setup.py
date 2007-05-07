#!/usr/bin/env python
"""
Only install lib as the tdl package
See setup_all.py to install modules as subpackages
"""

import distutils
from distutils.core import setup, Extension

from  lib.version import name,version,author,email,desc

package_dir = {'tdl': 'lib'}

#packages = ['tdl','tdl.TkPlotter','tdl.FileIO']
packages = ['tdl']

setup(
    name = name,
    version = version,
    author =  author,
    author_email = email,
    description = desc,
    package_dir = package_dir,
    packages = packages,
    license = 'Python',
    data_files  = [('bin',['tdl'])]
)
