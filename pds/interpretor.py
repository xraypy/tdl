#######################################################################
"""
T. Trainor (fftpt@uaf.edu), 10-2008  
Simple python interpretor

Modifications:
--------------

"""
#######################################################################
"""
Description:
Very simple set of classes that holds a SymbolTable
(Python namespace dictionary) and Executes python commands.

Note in general we should expect the following (as minimum)
methods/data in a SymbolTable and Interpretor

Interpretor::
.symbol_table          # 
.set_data(n,var)       #
.get_data(n)           #
.set_function(n,fcn)   #
.get_function(n)       #
.add_group(n)          #
.execute(s)            # run a single or list of statements
.evaluate(s)           # run a single statement, ret value
.get_input(p)          #

SymbolTable::
.add_symbol(n,val)     # 
.set_symbol(n,val)     #
.get_symbol(n)         #
.del_symbol(n)         #
.has_symbol(n)         #
.sym_type(n)           #
.list_symbols()        # return list of names

"""
###########################################################

import types
import sys
import time
try:
    import numpy as num
except:
    num = None

import util
import pyeval as eval

######################################################################
"""
Grouping of python types into our own categories
(see SymbolTable.sym_type):  
  'var' (DATATYPES)  - holds values or groupings
  'fnc' (FUNCTYPES)  - callable objects (functions and classes)
  'ins' (INSTTYPES)  - class instances 
  'mod' (MODTYPES)   - modules
  'oth' (OTHERTYPES) - some other python type we dont know how to group
"""
VARIABLE = "var"
FUNCTION = "fnc"
INSTANCE = "ins"
MODULE   = "mod"
OTHER    = "oth"

# General python 'variable' types
# i.e. things we normally consider
# as data or data containers
PYVARTYPES = [types.BooleanType,
              types.ComplexType,
              types.DictType,
              types.DictionaryType,
              types.FloatType,
              types.FileType,
              types.IntType,
              types.ListType,
              types.LongType,
              types.NoneType,
              types.StringType,
              types.StringTypes,
              types.TupleType,
              types.UnicodeType]

# Include numpy types
# NUMVARTYPES = [num.ArrayType] + num.typeDict.values() 
if num:
    NUMVARTYPES = [num.ndarray] + num.typeDict.values() 
    # VARTYPES are python + numpy 
    VARTYPES = PYVARTYPES + NUMVARTYPES 
else:
    VARTYPES = PYVARTYPES

# INSTTYPES are general user created objects
INSTTYPES = [types.InstanceType]

# Data includes normal variables, and instances
DATATYPES = VARTYPES + INSTTYPES

# Modules are modules
MODTYPES  = [types.ModuleType]

# GROUPTYPES lumps objects that
# contain generally user accessible
# attributes
GROUPTYPES = INSTTYPES + MODTYPES 

# FUNCTYPES includes functions and classes
# ie stuff you call as: x = func()
FUNCTYPES = [types.BuiltinFunctionType,
             types.BuiltinMethodType,
             types.ClassType,
             types.FunctionType,
             types.LambdaType,
             types.MethodType,
             types.TypeType,
             types.UnboundMethodType]

# Make sure these include numpy ufuncs
if num:
    FUNCTYPES = FUNCTYPES + [type(num.abs)]

"""
OTHERTYPES = [types.BufferType,
              types.CodeType,
              types.DictProxyType,
              types.EllipsisType,
              types.FrameType,
              types.GeneratorType,
              types.NotImplementedType,
              types.ObjectType,
              types.SliceType,
              types.TracebackType,
              types.TypeType,
              types.XRangeType,
              "<class 'site._Printer'>",
              "<class 'site._Helper'>"]
"""

######################################################################
######################################################################
class SymbolTable:

    def __init__(self):
        """
        Set up initial exectution name space.
        """
        self.data = eval.init_namespace()
        self._init_list()

    #########################################
    def add_symbol(self,name,value):
        #self.data.update({name:value})
        self.data['__tmp__'] = value
        s = "%s = __tmp__" % name
        eval.do_exec(s,self.data)
        del self.data['__tmp__']

    #########################################
    def set_symbol(self,name,value):
        #self.data[name]= value
        self.data['__tmp__'] = value
        s = "%s = __tmp__" % name
        eval.do_exec(s,self.data)
        del self.data['__tmp__']

    #########################################
    def del_symbol(self,name):
        #return self.data.pop(name)
        #del self.data[name]
        if self.has_symbol(name):
            s = "del %s" % name
            eval.do_exec(s,self.data)
        return

    #########################################
    """
    Note there are a number of ways to 'get' a symbol:
    - Assuming name is a key to the data dictionary.
            return self.name_space[name]
              or 
            return self.data.get(name)
    - Or assume x.y syntax and build using getattr.
            if name == None: return None
            if type(name) != types.StringType:
                raise "Error in get_symbol"
            if name.startswith('.'): return None
            top = None
            tmp = name.split('.')
            if len(tmp) > 0:
                top = self.data.get(tmp.pop(0))
            if top == None: return None
            for attr in tmp:
                top = getattr(top,attr)
            return top
    - Or build a str and execute.  This seems to be most general,
      ie can handle cases with arguments: data[1] etc..
            s = "try:\n"
            s = s+ "    __tmp__ = %s" % name
            s = s+ "except:"
            s = s +"    __tmp__ = None"
            eval.do_exec(s,self.data)
            tmp = self.data.get('__tmp__')
            return tmp
    """
    def get_symbol(self,name):
        s = "try:\n"
        s = s+ "    __tmp__ = %s\n" % name
        s = s+ "except:\n"
        s = s +"    __tmp__ = None\n"
        eval.do_exec(s,self.data)
        tmp = self.data.get('__tmp__')
        del self.data['__tmp__']
        return tmp

    #########################################
    def has_symbol(self,name):
        #if name in self.data.keys(): return True
        #else: return False
        if self.get_symbol(name) != None:
            return True
        else:
            return False

    #########################################
    def sym_type(self,name):
        """
        Return the type of symbol according to our grouping
        of actual python types.  
          var - holds values or groupings
          fnc - callable objects
          grp - class instances and modules
          oth - some other python type we dont know how to group
          unk - dont know what it is
        Return None if symbol doesnt exist
        """
        sym = self.get_symbol(name)
        if sym == None: return None

        ty = type(sym)
        if ty in VARTYPES:
            return VARIABLE
        elif ty in FUNCTYPES:
            if callable(sym):
                return FUNCTION
            else:
                return OTHER
        elif ty in INSTTYPES:
            return INSTANCE
        elif ty in MODTYPES:
            return MODULE
        else:
            return OTHER

    ######################################################
    def get_data_dict(self,name=None,_skip=True):
        """
        get toplevel data as a dictionary
        useful for saving data...
        """
        data = {}
        if name:
            d = self.get_symbol(name)
            if type(d) in DATATYPES:
                data[name] = d
        else:
            for key in self.data.keys():
                ignore = self._ignore(key,_skip=_skip)
                if not ignore:
                    ty = type(self.data[key])
                    if ty in DATATYPES:
                        data.update({key:self.data[key]})
        return data
    
    ######################################################
    def put_data_dict(self,data):
        """
        put top level data
        use to put saved data...
        """
        for key in data.keys():
            ty = type(data[key])
            if ty in DATATYPES:
                self.data.update({key:data[key]})
        return
    
    ######################################################
    def list_symbols(self,symbol=None,tunnel=True,_skip=True):
        self._init_list()
        self._list_symbol(symbol=symbol,tunnel=tunnel,_skip=_skip)
        self._list_data['var'].sort()
        self._list_data['fnc'].sort()
        self._list_data['ins'].sort()
        self._list_data['mod'].sort()
        self._list_data['oth'].sort()
        return self._list_data

    def list_symbols_all(self,symbol=None,tunnel=True,_skip=True):
        self._init_list()
        self._list_symbol(symbol=symbol,tunnel=tunnel,_skip=_skip)
        all = self._list_data['var']
        all = all + self._list_data['fnc']
        all = all + self._list_data['ins']
        all = all + self._list_data['mod']
        all = all + self._list_data['oth']
        all.sort()
        return all

    ######################################################
    def list_builtins(self,_skip=True):
        bins_list = {}
        bins_list['var']=[]
        bins_list['fnc']=[]
        bins_list['ins']=[]
        bins_list['mod']=[]
        bins_list['oth']=[]
        
        bins = self.get_symbol('__builtins__')
        if bins == None:
            raise "Error accessing __builtins__ in list_builtins"
        
        for name in bins.keys():
            ignore = self._ignore(name,_skip=_skip)
            if not ignore:
                ty = type(bins[name])
                if ty in VARTYPES:
                    bins_list['var'].append(name)
                elif ty in FUNCTYPES:
                    bins_list['fnc'].append(name)
                elif ty in INSTTYPES:
                    bins_list['ins'].append(name)
                elif ty in MODTYPES:
                    bins_list['mod'].append(name)
                else:
                    bins_list['oth'].append(name)
        bins_list['var'].sort()
        bins_list['fnc'].sort()
        bins_list['ins'].sort()
        bins_list['mod'].sort()
        bins_list['oth'].sort()
        return bins_list

    ######################################################
    def _init_list(self):
        """
        dictionary of lists:
        {'var':[variable names],
         'fcn':[function and method names],
         'ins':[object instance names],
         'mod':[module types]}
         'oth':[other py stuff]
        """
        self._list_data = {}
        self._list_data['var']=[]
        self._list_data['fnc']=[]
        self._list_data['ins']=[]
        self._list_data['mod']=[]
        self._list_data['oth']=[]
        
    def _list_symbol(self,symbol=None,tunnel=True,_skip=True):
        """
        list names of data and method atributes of a symbol
        if symbol is None we just pass to _list
        if symbol is not none we assume its the name of a "group"
        Note this will only "tunnel" on class instances
        """
        if symbol == None:
            self._list(tunnel=tunnel,_skip=_skip)
            return

        sym = self.get_symbol(symbol)
        if sym == None: return
        ty  = type(sym)

        if ty in GROUPTYPES:
            # Note: only tunnel on Instances.
            # tunneling on modules may result in
            # massive recursion....
            if ty not in INSTTYPES:
                tunnel = False
            lim = 0
            attrs = dir(sym)
            data = {}
            for attr in attrs:
                lim = lim + 1
                data.update({attr:getattr(sym,attr)})
                if lim > 1000:
                    print "Max attributes reached in _list_symbol..."
                    break
            self._list(data=data,pfx=symbol,tunnel=tunnel,_skip=_skip)
            return
        else:
            raise "Error in _list_symbol, passed argument is not a group"

    def _list(self,data=None,pfx=None,tunnel=True,_skip=True):
        if data==None:
            data = self.data
        lim = 0
        for key in data.keys():
            ignore = self._ignore(key,_skip=_skip)
            if not ignore:
                if pfx:
                    name = pfx + '.' + key
                else:
                    name = key
                ##
                ty = type(data[key])
                if ty in VARTYPES:
                    self._list_data['var'].append(name)
                elif ty in FUNCTYPES:
                    self._list_data['fnc'].append(name)
                elif ty in GROUPTYPES:
                    if ty in INSTTYPES:
                        self._list_data['ins'].append(name)
                    elif ty in MODTYPES:
                        self._list_data['mod'].append(name)
                    if tunnel:
                        self._list_symbol(symbol=name,tunnel=tunnel,_skip=_skip)
                        lim = lim + 1
                        if lim > 1000:
                            print "Max tunnel reached in _list..."
                            break
                else:
                    self._list_data['oth'].append(name)
        return
    
    def _ignore(self, key,_skip=True):
        if _skip == True:
            if key[0] == '_': ignore = True
            else: ignore = False
        else:
            ignore = False
        return ignore
    
##########################################################################
##########################################################################
USE_CODE = True
class Interpretor:
    
    def __init__(self,):

        self.symbol_table=SymbolTable()
        self.console     = None
        if USE_CODE:
            import code
            self.console = code.InteractiveConsole(self.symbol_table.data)
            
    ####################################################
    def execute(self,arg,**kws):
        """
        return 1 f command is completed
        return 2 if need more input
        """
        if USE_CODE:
            ret = self.console.push(arg)
            if ret == True:
                return 1
            else:
                return 0
        else:
            return eval.do_exec(arg,self.symbol_table.data)

    def execute_file(self,fname,**kws):
        return eval.do_execfile(fname,self.symbol_table.data)
        
    def evaluate(self,arg,**kws):
        return eval.do_eval(arg,self.symbol_table.data)

    ####################################################
    def set_data(self,name,val):
        if type(val) in DATATYPES:
            self.symbol_table.set_symbol(name,val)
        else:
            print "Error: %s not a valid data type (type=%s)" % (name,type(val))
    
    def get_data(self,name):
        val = self.symbol_table.get_symbol(name)
        if type(val) in DATATYPES:
            return val
        else:
            print "Error: %s not a valid data type (type=%s)" % (name,type(val))
        return None
    
    ####################################################
    def set_function(self,name,val):
        if type(val) in FUNCTYPES:
            if callable(val):
                self.symbol_table.set_symbol(name,val)
        else:
            print "Error: %s not a valid function type (type=%s)" % (name,type(val))
    
    def get_function(self,name):
        val = self.symbol_table.get_symbol(name)
        if type(val) in FUNCTYPES:
            return val
        else:
            print "Error: %s not a valid function type (type=%s)" % (name,type(val))
        return None
    
    ####################################################
    def add_group(self,name):
        val = Group()
        self.symbol_table.set_symbol(name,val)
    
##############################################
class Group:
    """generic group"""
    def __init__(self,):
        pass

##############################################################
##############################################################
if __name__ == "__main__":
    i = Interpretor()
    
    i.execute('x=10'); i.execute('print x')
    s = "import numpy as num\n"; i.execute(s)
    s = "xx = num.array([1,2,3])\n"; i.execute(s); i.execute('print xx')
    s = "from test.test import t\n"; i.execute(s)
    s = "tt = t()\nyy=tt.calc(10)\nprint yy"; i.execute(s)
    s = "qq = t()\nqq.ll=100\nqq.t=t()"; i.execute(s)
    
    ###
    print "\n----\n"
    d = i.symbol_table.listAll()
    print d

    print "\n----\n"
    d = i.symbol_table.listAll(group='qq')
    print d    
    
