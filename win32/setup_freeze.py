"""
setup.py script for cx_Freeze

Usage:
    python setup_freeze.py bdist_mac --iconfile=TDL.icns

then use disk utility to make a .dmg
"""

from cx_Freeze import setup, Executable

import os, sys
import wx
import wx.lib.agw.flatnotebook
import numpy, scipy, matplotlib
matplotlib.use('WXAgg')

from scipy.sparse.csgraph import _validation
import h5py

import Image

mpl_data_files = matplotlib.get_py2exe_datafiles()

resource_files = []
## also: use
# if getattr(sys, 'frozen', None) == 'macosx_app':
## in lib/site_config.py to specify plugin path

plugin_files = []
dll_files = []
# dll_files = [("larch/dlls/darwin/",
#               ['%s/dlls/darwin/%s' % (psrc, f) for f in os.listdir('%s/dlls/darwin' % psrc)]
#               )]

DATA_FILES = []
ICONFILE = 'TDL.icns'
import h5py
import h5py._objects, h5py._proxy, h5py.defs, h5py.utils
import tdl
from tdl.pds import *


pycard_incs = ['PythonCard', 'PythonCard.model', 'PythonCard.dialog']

import PythonCard.components
dirname = os.path.split(PythonCard.components.__file__)[0]
for comp in os.listdir(dirname):
    if comp.endswith('.py') and not comp.startswith('__'):
        pycard_incs.append('PythonCard.components.%s' % comp[:-3])

for tdl_mod in ('ana', 'ana_upgrade', 'geom', 'peak', 'specfile',
                 'spectra', 'sxrd', 'utils', 'xrd', 'xrf', 'xrr', 'xtab',
                 'xtal'):
    pycard_incs.append('tdl.modules.%s' % tdl_mod)


exe_opts = {'packages': ['wx', 'numpy', 'scipy', 'matplotlib'],
            'includes': ['Carbon', 'Carbon.Appearance', 'ConfigParser',
                         'Image', 'ctypes', 'fpformat', 'PythonCard',
                         'numpy', 'scipy', 'scipy.constants',
                         'scipy.fftpack', 'scipy.io.matlab.mio5_utils',
                         'scipy.io.matlab.streams', 'scipy.io.netcdf',
                         'scipy.optimize', 'scipy.signal',
                         'scipy.sparse.csgraph._validation', 'matplotlib',
                         'h5py', 'h5py._objects', 'h5py._proxy',
                         'h5py.defs', 'h5py.utils', 'tdl', 'tdl.pds',
                         'tdl.pds.pcgui', 'wx', 'wx._core', 'wx.lib',
                         'wx.lib.agw', 'wx.lib.agw.flatnotebook',
                         'wx.lib.agw.pycollapsiblepane',
                         'wx.lib.colourselect', 'wx.lib.masked',
                         'wx.lib.mixins', 'wx.lib.mixins.inspection',
                         'wx.lib.newevent', 'wx.py', 'wxversion', 'xdrlib',
                         'xml.etree', 'xml.etree.cElementTree'],
            'excludes': ['Tkinter', '_tkinter', 'Tkconstants', 'tcl',
                        '_imagingtk', 'PIL._imagingtk', 'ImageTk',
                        'PIL.ImageTk', 'FixTk''_gtkagg', '_tkagg',
                        'qt', 'PyQt4Gui', 'email', 'IPython'],
            #'iconfile': ICONFILE,
            }

exe_opts['includes'].extend(pycard_incs)
print pycard_incs

appname = 'PythonDataShell'
appvers = '0.1'
setup(name = appname,
      version = appvers,
      description = "GSECARS PythonDataShell",
      options = {'build_exe': exe_opts},
      data_files = mpl_data_files + dll_files, ##  + plugin_files
      executables = [Executable('runpds.py', base=None),
                     ])

print ' ===== '

contents = 'build/%s-%s.app/Contents' % (appname, appvers)
contents = contents.replace(' ', '\ ')

def sys(cmd):
    print ' >> ', cmd
    os.system(cmd)

sys("cp -pr TDL.icns  %s/." % contents)

try:
    os.makedirs("%s/Resources/tdl/" % contents)
except:
    pass
sys("cp -pr ../modules   %s/Resources/tdl/." % (contents))
#
# sys("cp -pr ../dlls/darwin/* %s/MacOS/." % contents)
