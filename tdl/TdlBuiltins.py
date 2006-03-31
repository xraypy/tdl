#!/usr/bin/env
# M. Newville Univ of Chicago (2005)
#
# --------------
# Modifications
# --------------
# 20-Mar-06 MN:  much rearranging, so that this modules holds all "builtins",
#             including some numeric stuff (sin,  pi, etc)
#         - eliminated 'delvar', as 'del' is now a real keyword
#         - renamed to TdlBuiltins (as this essentially holds all builtins"
#
# *1-28-06 T2:
#  moved numeric stuff to TdlNumLib
#
# * 1-12-06 T2:
#  - Renamed Lib -> TdlLib
#  - Reformat _func_ etc to match generic function loader in Evaluator
#  - Created a new tdl_import (import) to import/reload func modules
#  - Added a new tdl_delvar (delvar) function to delete variables and groups
#    works the same as setvar().  Note should del be added as a special keyword ?
#  - Add debug function to toggle tdl.debug
#  
##########################################################################
from Num import Num

import os
import sys
import types

from Util import show_list, show_more, datalen, unescape_string

title = "builtin library functions"

HelpBuiltins = """
  Builtin Functions: Overview

  additional help is available on these topics:
     help i/o      help on input/output (creating,writing,reading files)
     help os       help on interacting with the operating system (ls,path,...)
     help import   additional help on importing 
        
  Brief description of builtin functions:
     datagroup :  set default data group
         tdl> datagroup('mydata')

     funcgroup :  set default function group
         tdl> funcgroup('myfunctions')
         
     load      :  load and run a tdl program from external file
         tdl> load('program.tdl')
         
     import    :  import python modules
         tdl> import('module.py')        # imports (or re-imports) python module
         tdl> import()                   # re-imports all defined python modules
         
     eval      :  evaluate a tdl expression or statement
         tdl> eval('1/5')
         0.2
         tdl> eval(' x = sqrt(3)')
         tdl> print x
         1.73205080757 

     setvar    :  set a variable from a string, returns fully qualified name
         tdl> setvar('x',3)
         _main.x
         tdl> print x
         3.0

     type      :  returns the type of a variable

     read_ascii -
     debug      -
     abs        -
     max        -
     min        -
     len        -
     list       -

     strfind    -
     strsplit   -

     dictkeys   -
     dictvals   -
     dictitems  -

     cd         -
     pwd        -
     more       -
     path       -
     ls         -

"""


HelpIO = """
  Builtin Functions for Input/Output (creating,writing,reading files):

  = open(file,mode)
      tdl> f = open('myfile.txt','w')

  will open the file 'myfile.txt' for writing (and killing any existing file
  with that name!).  The value returned from open() is called the *file handle*
  and is a special value (a python object) that won't make sense as any other
  type.  But you can use a file handle to read, write, and close a file.

     tdl> write(f, 'a string\n')
     tdl> write(f, 'another line.\n')
     tdl> close(f)

  Note that writing a string to a file requires the newlines to be explicitly
  included.

  You can also read from an existing file if you open it in 'read mode':
     tdl> f = open('myfile.txt', 'r')
     tdl> x1 = readline(f)
     tdl> x2 = readline(f)
     tdl> print x1
     a string
     tdl> print x2
     another line

     tdl> close(f)
     
     
 
  
"""

HelpOS = """
  Builtin Functions for Interacting with the operating system

  ls

  more
  
  pwd

  cd

  !
  
"""

HelpImport = """
   Importing Python modules
"""

def _dictkeys(x):
    "return list of dictionary keys"
    if type(x) == types.DictType:
        return x.keys()
    else:
        return []

def _dictitems(x):
    "return list of dictionary items"
    if type(x) == types.DictType:
        return [list(i) for i in x.items()]
    else:
        return []
    
def _dictvals(x):
    "return list of dictionary keys"
    if type(x) == types.DictType:
        return x.values()
    else:
        return []    

def _len(x):
    "return length of data"
    return datalen(x)

def _list(x):
    "convert argument to a list type"
    try:
        return list(x)
    except TypeError:
        return x

def _strsplit(var,sep=' '):
    "split a string"
    if type(var) != types.StringType:
        print ' %s is not a string ' % var
        return None
    return var.split(sep)

def _strfind(var,sub):
    "string find"
    if type(var) != types.StringType:
        print ' %s is not a string ' % var
        return None
    return var.find(sub)

def _strchomp(var):
    "string find"
    # print 'var ', var, type(var), var[:-1]
    if type(var) != types.StringType:
        print ' %s is not a string ' % var
        return None
    return var[:-1]

def _strstrip(var,delim=None):
    "string find"
    if type(var) != types.StringType:
        print ' %s is not a string ' % var
        return None
    if delim==None:
        return var.strip()
    else: 
        return var.strip(delim)

def _ls(arg = '*'):
    " return list of files in the current directory "
    from glob import glob
    if arg != None and type(arg)==types.StringType:
        return glob(arg)
    else:
        return os.listdir(os.getcwd())

def _ls_cmdout(x,ncol=None):
    " output for ls "
    return show_list(x,ncol=ncol)

def _cwd(x=None):
    "return current working directory"
    return os.getcwd()

def _cd(name):
    "change directorty"
    os.chdir(name)
    return os.getcwd()    

def _more(name,pagesize=24):
    "list file contents"
    try:
        f = open(name)
        l = f.readlines()
        f.close()
        show_more(l,filename=name,pagesize=pagesize)
    except IOError:
        print "cannot open file: %s." % name
        return
    

def _type(x):
    "return data type of data"
    t = type(x)
    # print 'XX TDL TYPE ', x, type(x)
    typecodes = {types.StringType:'string',
                 types.IntType: 'int',
                 types.LongType:'int',
                 types.FloatType:'float',
                 types.ComplexType:'complex',
                 types.ListType:'list',
                 types.DictType:'dict',
                 Num.ArrayType:'array'}
    
    if t in typecodes.keys(): return typecodes[t]
    return 'unknown'             


def _path(add=None):
    "modify python path"
    if add == None:
        for l in sys.path:
            print l
    else:
        sys.path.append(add)
    return

def tdl_open(filename,mode='r',tdl=None,**kw):
    " open a file "
    return open(filename,mode=mode)

def tdl_close(file,tdl=None,**kw):
    " close a file "
    if type(file) == types.FileType:
        return file.close()
    else:
        return False


def tdl_write(file,s,tdl=None,**kw):
    " write a string to a file "
    if type(file) == types.FileType:        
        return file.write(unescape_string(s))

def tdl_flush(file,tdl=None,**kw):
    " flush a file "
    if type(file) == types.FileType: return file.flush()

def tdl_read(file,size=None,tdl=None,**kw):
    " read from a file "
    if type(file) == types.FileType: return file.read(size=size)

def tdl_readline(file,size=None,tdl=None,**kw):
    " readline a file "
    if type(file) == types.FileType: return file.readline(size=size)    

def tdl_readlines(file,tdl=None,**kw):
    " read all lines of text from a file "
    if type(file) == types.FileType: return file.readlines()

def tdl_seek(file,offset,whence=None,tdl=None,**kw):
    " read all lines of text from a file "
    if type(file) == types.FileType: return file.seek(offset,whence=whence)        

def tdl_tell(file,offset,whence=None,tdl=None,**kw):
    " read all lines of text from a file "
    if type(file) == types.FileType: return file.tell()

def tdl_set_debug(debug=None,tdl=None,**kw):
    if tdl == None:
        return None
    if debug == None:
        debug = not tdl.debug
        
    tdl.set_debug(debug)
    return None

def tdl_load(fname, tdl=None,debug=False,**kw):
    " load file of tdl code"
    if tdl == None:
        if debug: print 'cannot run file %s ' % fname
        return None
    if debug: print 'load .... ', fname
    if not os.path.exists(fname):
        print 'file error: cannot find file %s ' % fname
        return None        
    tdl.load_file(fname)
    tdl.run()
    if debug: print 'load done.'

def tdl_import(lib='', tdl=None,debug=False,reloadAll=False,clearAll=False,**kw):
    """
    import python modules that define tdl functions,
    import()               # re-imports all previously defined modules
    load('x.py')           # re-imports all modules and include new module x.py
    import(clearAll=True)  # re-imports modules AND clears all data  
    """
    if tdl == None:
        if debug: print 'cannot load function modules ' 
        return None
    if debug: print 'loading function modules.... '
    if lib=='' or reloadAll:
        tdl.symbolTable.initialize(clearAll=clearAll)
        
    tdl.symbolTable.import_lib(lib)
    if debug: print 'import done.'

def tdl_eval(expr, tdl=None,debug=False,**kw):
    " evaluate tdl expression"
    if tdl == None:
        if debug: print 'cannot eval %s ' % expr
        return None
    return tdl.do_eval(expr)

def tdl_setvar(name,val,tdl=None,group=None,debug=False,**kws):
    "set default group"
    # print 'This is tdl setvar ', name, val, tdl, group, kws
    if tdl == None:
        if debug: print 'cannot setvar %s ' % expr
        return None
    name.strip()
    #idot = name.find('.')
    #if idot > 1:
    #    group = name[:idot]
    #    name  = name[idot+1:]
    
    return '.'.join(tdl.symbolTable.setVariable(name,val))

def tdl_delvar(name,tdl=None,group=None,debug=False,**kw):
    "delete a variable or group"
    if tdl == None:
        if debug: print 'cannot setvar %s ' % expr
        return None
    name.strip()
    xx = name.split('.')
    if len(xx) == 1:
        return tdl.symbolTable.deleteSymbol(name,group=group)
    elif len(xx) == 2 and xx[1] == '':
        return tdl.symbolTable.deleteGroup(xx[0])
    else:
        return tdl.symbolTable.deleteSymbol(xx[1],group=xx[0])

def tdl_group2dict(gname=None,tdl=None,debug=False,**kw):
    "convert all data in a group to a single dictionary"
    if tdl == None:
        if debug: print 'cannot setgroup %s ' % expr
        return None
    if gname == None: gname = tdl.symbolTable.getDataGroup()
    gname = gname.strip()
    dict = {}
    dat = tdl.symbolTable.getAllData(group=gname)
    for i in dat:
        if i.type not in ('pyfunc','defpro'):
            dict[i.name] = i.value
    return dict


def tdl_setdatagroup(gname=None,tdl=None,debug=False,**kw):
    "set default group"
    if tdl == None:
        if debug: print 'cannot setgroup %s ' % expr
        return None
    if gname == None:
        return tdl.symbolTable.getDataGroup()
    g = gname.strip()
    return tdl.symbolTable.setDataGroup(g)

def tdl_setfuncgroup(gname=None,tdl=None,debug=False,**kw):
    "set default group"
    if tdl == None:
        if debug: print 'cannot setgroup %s ' % expr
        return None
    if gname == None:
        return tdl.symbolTable.getFuncGroup()
    g = gname.strip()
    return tdl.symbolTable.setFuncGroup(g)

def tdl_func_as_cmd(name,tdl=None):
    "allow functions to act as commands"
    if tdl == None:
        if debug: print 'cannot setgroup %s ' % expr
        return None
    if tdl.symbolTable.hasFunc(name):
        sym = tdl.symbolTable.getSymbol(name)
        sym.as_cmd = True
    return


def tdl_read_ascii(fname, group=None, tdl=None,debug=False, **kw):
    " read ascii file of tdl code"
    from ASCIIFile import ASCIIFile
    if tdl == None:
        if debug: print 'cannot read file %s ' % fname
        return None
    if debug: print 'reading.... ', fname
    if not os.path.exists(fname):
        print 'read_ascii: cannot find file %s ' % fname
        return None        
    try:
        f = ASCIIFile(fname)
    except:
        print 'read_ascii: error reading file %s ' % fname        
        return None
    
    # save current group name
    savegroup = tdl.symbolTable.getDataGroup()
    if group == None:
        group = savegroup
    else:
        group = tdl.symbolTable.setDataGroup(group)
    tdl.symbolTable.setDataGroup(group)
    
    tdl.symbolTable.setVariable('titles', f.get_titles())
    tdl.symbolTable.setVariable('column_labels', f.labels)        
    ret = [group]
    for i in f.labels:
        ret.append(i)
        tdl.symbolTable.setVariable(i, f.get_array(i))
    if debug: print 'read done.'
    # return default group to original
    tdl.symbolTable.setDataGroup(savegroup)
    return ret
    

def _int(x):
    "wrap builtin int or Numpy astype(int)"
    if type(x) == Num.ArrayType:
        return x.astype(int)
    else:
        return int(x)

def _float(x):
    "wrap builtin float or numpy astype(float)"
    if type(x) == Num.ArrayType:
        return x.astype(float)
    else:
        return float(x)

def _complex(x):
    "wrap builtin complex or numpy astype(complex)"
    if type(x) == Num.ArrayType:
        return x.astype(complex)
    else:
        return complex(x)


def _range(x,stop=None,step=None,shape=None,dtype='d'):
    """create an array of evenly spaced values:
        range(x, [stop=stop, [step=step, [shape=shape, [dtype='d']]]])

    Thus,
       range(10)
    returns [0. 1. 2. 3. 4. 5. 6. 7. 8. 9.] (the first 10 numers, as double precision)
    other variations: 
       range(2,5)               -> [2. 3. 4.]
       range(2,10,2)            -> [2. 4. 6. 8.]
       range(2,10,2,dtype='i')  -> [2 4 6 8]

       range(20,shape=[5,4])    ->  [[  0.   1.   2.   3.]
                                     [  4.   5.   6.   7.]
                                     [  8.   9.  10.  11.]
                                     [ 12.  13.  14.  15.]
                                     [ 16.  17.  18.  19.]]
    """
    if stop==None and step==None:         t = Num.arange(x,dtype=dtype)
    elif stop != None and step == None:   t = Num.arange(x,stop,dtype=dtype)
    elif stop != None and step != None:   t = Num.arange(x,stop,step,dtype=dtype)
    elif stop == None and step != None:   t = Num.arange(0, x,step,dtype=dtype)

    if shape is None: return t
    reshaped = False
    if type(shape) in  (types.ListType,types.TupleType):
        npts = 1
        for i in shape: npts = npts*int(i)
        if npts == len(t):
            try:
                t = Num.array(t)
                t.shape = tuple(shape)
                return t
            except ValueError:
                pass
    # if we got here, there's a value error
    raise ValueError, ' could not re-shape array to shape = %s ' % repr(shape)
            

#################################################################
# Load the functions
#################################################################

# constants to go into _builtin name space
_consts_ = {"True": True, "False":False, "None":None,
            "pi":Num.pi, "e":Num.e}            

_help_  = {'builtins': HelpBuiltins,
           'i/o':HelpIO,
           'os': HelpOS,
           'import': HelpImport,
           }

#'_builtin.del':(tdl_delvar,None),

# functions to add to namespace
_func_ = {'_builtin.load':(tdl_load, None),
          '_builtin.import':(tdl_import, None),
          '_builtin.eval':(tdl_eval,None),
          '_builtin.setvar':(tdl_setvar,None),
          '_builtin.datagroup':(tdl_setdatagroup, None),
          '_builtin.funcgroup':(tdl_setfuncgroup, None),
          '_builtin.ascmd':(tdl_func_as_cmd,None),
          '_builtin.read_ascii':(tdl_read_ascii,None),
          '_builtin.debug':(tdl_set_debug,None),
          "_builtin.cd":(_cd,None),
          "_builtin.pwd":(_cwd,None),
          "_builtin.more":(_more,None),
          "_builtin.path":(_path,None),
          "_builtin.ls":(_ls,_ls_cmdout),
          "_builtin.abs":(abs,None),
          "_builtin.max":(max,None),
          "_builtin.min":(min,None),
          "_builtin.len":(_len,None),
          "_builtin.list":(_list,None),
          "_builtin.strfind":(_strfind,None),
          "_builtin.strsplit":(_strsplit, None),
          "_builtin.strstrip":(_strstrip, None),
          "_builtin.strchomp":(_strchomp, None),
          "_builtin.type":(_type,None),
          "_builtin.dictkeys":(_dictkeys,None),
          "_builtin.dictvals":(_dictvals,None),
          "_builtin.dictitems":(_dictitems,None),
          "_builtin.array":Num.array,
          "_builtin.sin":Num.sin,
          "_builtin.cos":Num.cos,
          "_builtin.tan":Num.tan,          
          "_builtin.exp":Num.exp,
          "_builtin.ln":Num.log,
          "_builtin.log":Num.log,
          "_builtin.log10":Num.log10,
          "_builtin.sqrt":Num.sqrt,
          "_builtin.int": _int,
          "_builtin.float": _float,
          "_builtin.complex": _complex,
          "_builtin.range":_range,
          "_builtin.arange":_range,
          "_builtin.open":tdl_open,
          "_builtin.close":tdl_close,
          "_builtin.write":tdl_write,
          "_builtin.read":tdl_read,
          "_builtin.readlines":tdl_readlines,
          "_builtin.group2dict":tdl_group2dict,
          
          }


if __name__ == '__main__':
    print help
