#!/usr/bin/env python
"""
Setup file to install tdl into site packages.
  >>python setup.py install
Or to prevent copying to site-packages
  >>python setup.py build

When installed the resulting package should look like:
    site-packages/tdl/doc
    site-packages/tdl/lib
    site-packages/tdl/modules
    site-packages/tdl/pds
    site-packages/tdl/scripts

"""
#########################################################

import distutils
from distutils.core import setup, Extension

#from version import name,version,author,email,desc
version = '0.4'
name    = 'tdl'
author  = "pyxrd"
email   = "xxx"
desc    = "tdl"

#########################################################
package_dir  = {'tdl':'',
                'tdl.doc':'doc',
                'tdl.lib':'lib',
                'tdl.modules':'modules',
                'tdl.pds':'pds',
                'tdl.scripts':'scripts'}
packages     = ['tdl',
                'tdl.doc',
                'tdl.lib',
                'tdl.modules',
                'tdl.pds',
                'tdl.scripts']

# misc package data
package_data = {}
package_data.update({'tdl.doc':['LICENSE.txt']})
package_data.update({'tdl.scripts':['pds']})
package_data['tdl.scripts'].append('pds.bat')
package_data['tdl.scripts'].append('spectohdf.py')
                
#pds
packages.append('tdl.pds.pcgui')
packages.append('tdl.pds.menu')
package_data.update({'tdl.pds':['startup.pds']})

# tdl.lib
package_data.update({'tdl.lib':['_xrr.dll']})
package_data['tdl.lib'].append('_hello.lib')
package_data['tdl.lib'].append('gsl.dll')
package_data['tdl.lib'].append('gsl.lib')
package_data['tdl.lib'].append('gslcblas.dll')

# modules
packages.append('tdl.modules.ana_upgrade')
packages.append('tdl.modules.ana')
packages.append('tdl.modules.geom')
packages.append('tdl.modules.peak')
packages.append('tdl.modules.specfile')
packages.append('tdl.modules.spectra')
packages.append('tdl.modules.sxrd')
packages.append('tdl.modules.utils')
packages.append('tdl.modules.utils.mpfit')
packages.append('tdl.modules.utils.files')
packages.append('tdl.modules.xrd')
packages.append('tdl.modules.xrf')
packages.append('tdl.modules.xrr')
packages.append('tdl.modules.xtab')
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
    data_files  = []
)
