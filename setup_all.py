#!/usr/bin/env python
"""
This setup file defines the directory lib as tdl
and the directory modules as tdl.modules.

The resulting installed root package name is tdl
and include all the modules defined in lib (ie lib/*.py) 
We also define the root package to include a sub package tdl.modules
which inlcudes all the modules defined in modules (ie modules/*.py).

Additional packages under modules/ must be given 
explicitly - see append below.
"""

import distutils
from distutils.core import setup, Extension

# Import tdl version information
from  lib.version import name,version,author,email,desc


### Define lib and modules as tdl and tdl.modules respectivley
package_dir = {'tdl': 'lib','tdl.modules': 'modules'}
packages = ['tdl','tdl.modules']
package_data = {'tdl.modules':['README.modules','startup.tdl']}

### Add-on modules
packages.append('tdl.modules.TkPlotter')
packages.append('tdl.modules.ASCIIFile')
packages.append('tdl.modules.wxGUI')
packages.append('tdl.modules.xray')
packages.append('tdl.modules.xray.CarsMca')
packages.append('tdl.modules.xray.Spec')
packages.append('tdl.modules.xray.ScanData')
package_data.update({'tdl.modules.xray':['spec_scripts.tdl','xrf_scripts.tdl']})

### Additional data files 
data_files = [('bin',['tdl'])]

### call the setup command
setup(
    name = name,
    version = version,
    author =  author,
    author_email = email,
    description = desc,
    package_dir = package_dir,
    packages = packages,
    package_data = package_data,
    license = 'Python',
    data_files  = data_files
)
