#!/usr/bin/env python

import distutils
from distutils.core import setup, Extension

# Import tdl version information
from  lib.version import name,version,author,email,desc

# main package directory under .../python/.../site-packages will be named tdl
package_dir = {'tdl': ''}

# the packages to be installed include tdl (py files in same dir as this script)
# and the lib and modules sub directories
packages = ['tdl','tdl.lib','tdl.modules']
package_data = {'tdl.modules':['README.modules','startup.tdl']}

### Add-on modules
packages.append('tdl.modules.TkPlotter')
packages.append('tdl.modules.ASCIIFile')
packages.append('tdl.modules.GUI')
packages.append('tdl.modules.xray')
packages.append('tdl.modules.xray.CarsMca')
packages.append('tdl.modules.xray.Spec')
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
