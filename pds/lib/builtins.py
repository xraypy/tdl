"""
Builtin functions for pds

Authors/Modifications:
----------------------
- T. Trainor, M. Newville
- Modified from tdl r259 for use with pds shell program

"""
#####################################################################

import os
import sys
import types
import time

from shellutil import show_list, show_more, datalen
from shellutil import set_path, unescape_string, list2array
from shellutil import mod_import
from shellutil import Group

#####################################################################
class PdsBuiltins:
    """
    Builtin methods / fcns
    """
    #################################################################
    def __init__(self,):
        """ init """
        self.module_list = []

    #################################################################
    def group(self,**kw):
        """
        Create a new group

        Examples:
        -------- 
        >>g = group()
        >>g = group(x=x,y=data['y'])
        In the second example 'g.x' and 'g.y' will
        be assigned from the key word arguments
        Once you have a group you can modify it on the fly
        >>g.zz = 100
        >>g.x = num.sin(10)
        """
        grp = Group()
        if len(kw) == 0:
            return grp
        for key in kw.keys():
            setattr(grp,key,kw[key])
        return grp

    #################################################################
    def ls(self,arg= '.',ncol=1,width=72):
        """
        List directory contents

        Parameters:
        ----------
        * arg is the path
        * ncol is the the number of columns to display
        * textwidth is the width for display
        """
        ret = self._ls(arg=arg)
        print show_list(ret,ncol=ncol,textwidth=width)
        print ""

    #################################################################
    def _ls(self,arg= '.'):
        """ Return list of files in the given directory """
        from glob import glob
        arg.strip()
        if type(arg) != types.StringType or len(arg)==0: arg = '.'
        if os.path.isdir(arg):
            ret = os.listdir(arg)
        else:
            ret = glob(arg)
        if sys.platform == 'win32':
            for j in range(len(ret)):
                ret[j] = ret[j].replace('\\','/')
        return ret

    #################################################################
    def pwd(self,):
        """Print current working directory"""
        print self._cwd()
        print ""

    #################################################################
    def _cwd(self,x=None):
        """Return current working directory"""
        ret = os.getcwd()
        if sys.platform == 'win32':
            ret = ret.replace('\\','/')
        return ret

    #################################################################
    def cd(self,name=None):
        """Change current working directory"""
        self._cd(name)
        
    #################################################################
    def _cd(self,name=None):
        """Change directory"""
        if name == None:
            self.pwd()
            return  
        if type(name) == types.StringType:
            name = name.strip()
        else:
            if hasattr(name,'__name__'):
                name = name.__name__
            else:
                name = str(name)
        if name:
            try:
                os.chdir(name)
            except:
                print "Directory '%s' not found" % name
        ret = os.getcwd()
        if sys.platform == 'win32':
            ret = ret.replace('\\','/')
        return ret

    #################################################################
    def more(self,name,pagesize=24):
        """
        List file contents
        """
        try:
            f = open(name)
            l = f.readlines()
            f.close()
            show_more(l,filename=name,pagesize=pagesize)
        except IOError:
            print "cannot open file: %s." % name
            return

    #################################################################
    def path(self,pth=None,recurse=False,**kw):
        """
        Modify or show python path

        * if called with no arguments show sys.path
        * if passed the pth argument add it to sys.path
          (assume pth is a string path name)
        * if recurse=True all subdirectories (with names
           starting with an alpha char) are added
        """
        ret = self._path(pth=pth,recurse=recurse,**kw)
        if ret: print ret
        
    #################################################################
    def _path(self,pth=None,recurse=False,**kw):
        """
        Modify or show python path
        """
        if not pth:
            return show_list(sys.path)
        else:
            set_path(pth=pth,recurse=recurse)
        return None

    #################################################################
    def rimport(self,module=None,**kw):
        """
        Import/reimport python modules

        Examples:
        --------
        >>x = rimport('x')
        Imports the new module 'x.py' similiar to a regular python import
        however, we keep track of this module in list for reimports
        
        >>rimport()     # re-import all rimported modules
        >>rimport(x)    # just re-import x
        >>rimport('x')  # just re-import x
        
        This re-imports previously defined modules so
        any recent changes are now loaded
        """
        def _import(module):
            mod = mod_import(module)
            if type(mod) == types.ModuleType:
                if mod not in self.module_list: 
                    self.module_list.append(mod)
            return mod

        if module:
            return _import(module)
        else:
            for m in self.module_list:
                mod = _import(m)
                
    #################################################################
    def source(self,object,**kw):
        """
        Show the source code of a python object

        Examples:
        ---------
        >>source(my_function)  
        """
        import inspect
        src = inspect.getsource(object)
        for l in src.split('\n'): print l

    #################################################################
    def info(self,item):
        """
        Print useful information about item.

        Examples:
        ---------
        >>info(my_function)  
        """
        if hasattr(item, '__name__'):
            print "NAME:    ", item.__name__
        if hasattr(item, '__class__'):
            print "CLASS:   ", item.__class__.__name__
        print "ID:      ", id(item)
        print "TYPE:    ", type(item)
        print "VALUE:   ", repr(item)
        print "CALLABLE:",
        if callable(item):
            print "Yes"
        else:
            print "No"
        if hasattr(item, '__doc__'):
            doc = getattr(item, '__doc__')
        doc = doc.strip()   
        firstline = doc.split('\n')[0]
        print "DOC:     ", firstline
        
#####################################################################
# Load the functions on import
#####################################################################
__pdsbuiltins__ = PdsBuiltins()

