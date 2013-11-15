import os
import sys
from platform import uname, architecture

TDL_TOPDIR = '.'

def add2path(envvar='PATH', dirname='.'):
    """add specified dir to begninng of PATH and
    DYLD_LIBRARY_PATH, LD_LIBRARY_PATH environmental variables,
    returns previous definition of PATH, for restoration"""
    sep = ':'
    if os.name == 'nt':
        sep = ';'
    oldpath = os.environ.get(envvar, '')
    if oldpath == '':
        os.environ[envvar] = dirname
    else:
        paths = oldpath.split(sep)
        paths.insert(0, os.path.abspath(dirname))
        os.environ[envvar] = sep.join(paths)
    return oldpath

def get_dlldir():
    system, node, release, version, mach, processor = uname()
    arch = architecture()[0]
    dlldir = None
    suff = '32'
    if arch.startswith('64'):  suff = '64'
    if os.name == 'nt':
        return 'win%s' % suff
    elif system.lower().startswith('linux'):
        return 'linux%s' % suff
    elif system.lower().startswith('darwin'):
        return 'darwin'
    return ''

def get_dll(libname):
    """find and load a shared library"""
    _paths = {'PATH': '', 'LD_LIBRARY_PATH': '', 'DYLD_LIBRARY_PATH':''}
    _dylib_formats = {'win32': '%s.dll', 'linux2': 'lib%s.so',
                      'darwin': 'lib%s.dylib'}
    thisdir = os.path.join(TDL_TOPDIR, 'dlls', get_dlldir())
    thisdir = os.path.abspath(thisdir)
    dirs = [thisdir]

    loaddll = ctypes.cdll.LoadLibrary

    if sys.platform == 'win32':
        loaddll = ctypes.windll.LoadLibrary
        dirs.append(TDL_TOPDIR)

    if hasattr(sys, 'frozen'): # frozen with py2exe!!
        dirs.append(os.path.dirname(sys.executable))

    for key in _paths:
        for d in dirs:
            _paths[key] = add2path(key, d)

    # normally, we expect the dll to be here in the dlls tree
    # if we find it there, use that one
    fname = _dylib_formats[sys.platform] % libname
    dllpath = os.path.join(thisdir, fname)
    if os.path.exists(dllpath):
        return loaddll(dllpath)

    # if not found in the dlls tree, try your best!
    return loaddll(ctypes.util.find_library(libname))

