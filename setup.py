#!/usr/bin/env python
"""
Setup file to install tdl into site packages.
  >>python setup.py install
Or to prevent copying to site-packages
  >>python setup.py build

When installed the resulting package should look like:
  site-packages/tdl/astlib
  site-packages/tdl/pds
  site-packages/tdl/modules

"""
#########################################################

import distutils
from distutils.core import setup, Extension

#from  lib.version import name,version,author,email,desc
version = '0.9'
name    = 'tdl'
author  = "Newville and Trainor"
email   = "xxx"
desc    = "tdl"

#########################################################
package_dir  = {'tdl.astlib':'astlib',
                'tdl.pds':'pds',
                'tdl.modules':'modules', }
packages     = ['tdl.astlib','tdl.pds','tdl.modules']
package_data = {'tdl.pds':['startup.pds']}

### Add-on modules
#wrap and lib
packages.append('tdl.modules.wrap')
package_data.update({'tdl.modules.wrap':['libs/_ref.dll']})
package_data['tdl.modules.wrap'].append('libs/core.lib')
package_data['tdl.modules.wrap'].append('libs/gsl.dll')
package_data['tdl.modules.wrap'].append('libs/gsl.lib')
package_data['tdl.modules.wrap'].append('libs/gslcblas.dll')
#mpcutils
packages.append('tdl.modules.mpcutils')
packages.append('tdl.modules.mpcutils.mpfit')
#fmtfiles
packages.append('tdl.modules.fmtfiles')
#wxgui
packages.append('tdl.modules.wxgui')
#xray
packages.append('tdl.modules.xray')
packages.append('tdl.modules.xray.detector')
packages.append('tdl.modules.xray.scandata')
packages.append('tdl.modules.xray.xrd')
packages.append('tdl.modules.xray.xrf')
packages.append('tdl.modules.xray.xrr')
packages.append('tdl.modules.xray.xtab')
#xtal
packages.append('tdl.modules.xtal')

### call the setup command
setup(
    name = name,
    version = version,
    author =  author,
    author_email = email,
    description = desc,
    license = 'Python',
    package_dir = package_dir,
    packages = packages,
    package_data = package_data,
    data_files  = [('bin',['tdl'])]
)
