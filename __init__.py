"""
Import/setup for the tdl modules

We put all the module subdirectories on the path.
Therefore, if you do the following (assuming the tdl
directory is somewhere on your PYTHONPATH):

>>import tdl

Then all the directories within the tdl/modules directory
will be added to your path.  Therefore, you can import
from any submodule without having to use the tdl.modules.whatever
syntax.  in other words the following should work.

>>import whatever

"""
try:
    from pds.lib.shellutil import set_path
    import os, sys
    p = os.path.dirname(__file__)
    set_path(p,recurse=False)
    p = os.path.join(p,'modules')
    set_path(p,recurse=True)
    del set_path, p
except:
    print "Warning: Error setting tdl path"
    pass
