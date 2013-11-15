# main exe
import os, sys
import wx
import numpy, scipy, h5py
from tdl.pds import shell

if hasattr(sys, 'frozen'):
    if os.name == 'nt':
        try:
            tdir, exe = os.path.split(sys.executable)
            toplevel, bindir = os.path.split(tdir)
            sys.path.append(os.path.abspath(tdir))
            sys.path.append(os.path.join(tdir, 'modules'))
        except:
            pass
    elif sys.platform.lower().startswith('darwin'):
        tdir, exe = os.path.split(sys.executable)
        toplevel, bindir = os.path.split(tdir)
        sys.path.append(os.path.join(toplevel, 'Resources', 'tdl'))
        sys.path.append(os.path.join(toplevel, 'Resources', 'tdl', 'modules'))
        sys.path.append(os.path.join(toplevel, 'Resources', 'tdl', 'pds'))

shell.main(use_wx=True)

