#!/usr/bin/env python

import distutils
from distutils.core import setup, Extension

from  lib.version import name,version,author,email,desc

# still not sure where to look for
# domain-specific modules and user modules
# environmental variables??

package_dir = {'tdl': 'lib'}

## , 'tdl.modules': 'modules'}

packages = ['tdl','tdl.TkPlotter','tdl.FileIO']
## ,'tdl.modules','tdl.modules.GUI','tdl.modules.IO']


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
