# M. Newville Univ of Chicago (2006)
#########################################################################

__version__  = '0.7'

# symbol Groups that must always be present, and starting search order
initGroups    = ['_sys','_math','_plot','_builtin']
init_loadLibs = ['TdlBuiltins','TdlNumLib']
_TopGroupName = '_'

import types
import os
import sys
import re
import random
import inspect
import exceptions
from copy import deepcopy

import Num

from Util import set_path, list2array
from Util import  show_more
from string import ascii_lowercase, digits

isValidName = re.compile(r'^[a-zA-Z_\$&@][a-zA-Z_\$&@0-9]*$').match

def randomName(n=8):
    return ''.join([(ascii_lowercase + digits)[random.randrange(36)] for i in range(n)])

class SymbolError(exceptions.Exception):
    def __init__(self,error,descr = None,node = None):
        self.error = error
        self.descr = descr
    def __repr__(self):
        return "%s: %s" % (self.error, self.descr)
    __str__ = __repr__

def datalen(x):
    "return length of data for many datatypes"
    try:
        return len(x)
    except:
        return 1

def set_path(p):
    print 'Dummy Set Path ', p


def list2array(x):
    """ attempt to convert a list to a NumPy array.
    Returns original list if the conversion is not possible"""

    def _IsNumeric(x):
        """returns whether a sequence of value (potentially nested list)
        can be coerced to numpy array of all numbers.
        """
        for i in x:
            if not (isinstance(i, (int,long,float,complex)) or (isinstance(i,(list,tuple)) and _IsNumeric(i))):
                return False
        return True

    if isinstance(x,(list,tuple)):
        try:
            if _IsNumeric(x): return Num.array(x)
        except ValueError:  # may be a list of variable-length lists or something
            pass 
    return x
    
def splitName(name,check=True):
    parts = name.split('.')
    if check:
        for i in parts:
            if not isValidName(i):
                msg = 'invalid name "%s"' % (name)
                if len(n)>1: msg = 'invalid name "%s" (error at "%s")' % (name,i)
                raise SymbolError, msg
    return parts
    
class symTypes:
    variable  = 'variable'
    group     = 'group'
    pyobj     = 'python object'
    file      = 'file'
    defvar    = 'defined variable'
    pyfunc    = 'python function'
    defpro    = 'tdl procedure'
    Data      = (variable,group,pyobj,file,defvar)
    Funcs     = (pyfunc,defpro)
    All       = Data + Funcs

def isSymbol(obj): return isinstance(obj, Symbol)
def isGroup(obj):  return isinstance(obj, Group)

class Group(dict):
    STATCODES = ('normal','frozen','nodelete','delete')
    def __init__(self,name='',filename=None, status='normal'):
        dict.__init__(self)
        self.name = name
        self.filename = filename
        if status not in self.STATCODES: status = 'normal'
        self.status = status
        
    def _symbols(self):
        s = []
        for k,v in self.items():
            if isSymbol(v): s.append(k)
        return s

    def _hasSymbol(self,sym):
        return (self.has_key(sym) and isSymbol(self[sym]))

    def _groups(self):
        s = []
        for k,v in self.items():
            if isGroup(v): s.append(k)
        return s

    def _hasGroup(self,sym):
        return (self.has_key(sym) and isGroup(self[sym]))

class Symbol:
    """
    Basic container for all tdl symbols. 
    Each symbol has the following components:
 
       name       Simple Name used to access symbol
       value      Python object for value
       type       Tdl variable type (not python type! a member of symTypes above)
       constant   boolean for whether value is fixed.
       desc       optional description (stores text formula for defvar)
       help       help string (say, for defpro)
       code       Tdl code for defpro
       args       args for defpro
       kws        keywords for defpro
       module     module namespace for defpro
       cmd_out    method for processing output when a pyfunc of 
                  defpro is run as a command.
    Notes:
       a Simple Name matches [a-zA-Z_\$][a-zA-Z_\$0-9]*
       Of course, a Symbol may be part of a hiearchy of Groups,
       and so may have a much longer Full Name.
    """
    name     = ''
    constant = False
    def __init__(self, name, value=None, constant=False,
                 type=None, desc=None, help=None, cmd_out=None,
                 code=None, args=None, kws=None, module=None):
        self.constant = False
        self.name     = name
        self.cmd_out  = cmd_out
        self.code     = code
        self.desc     = desc or ''
        self.help     = help or ''
        self.args     = args or []
        self.kws      = kws  or {}
        self.module   = module

        if type is None: type = symTypes.variable
        
        # check that type is valid
        if type not in symTypes.All:
            raise SymbolError, 'cannot add symbol "%s" with type=%s' % (name,type)

        self.type = type
        self.constant = constant
        # note that value is set last, as this may have side effects
        self.value    = value


    def getHelp(self):
        return self.help

    def getCode(self):
        return deepcopy(self.code)

    def __call__(self, *args,**kws):
        if self.type != symTypes.pyfunc: return self.value

        x = self.kws.copy() ; x.update(kws)
        return self.value(*args,**x)

    def cmdout(self,val,**kws):
        if val is None: return None
        if self.cmd_out is None: return str(val)

        x = self.kws.copy() ; x.update(kws)
        return self.cmd_out(val,**x)
         
    def __setattr__(self, attr, val):
        """ here to prevent re-setting of constants"""
        if self.constant and attr == 'value':
            raise SymbolError,  'cannot set value of constant %s' % (self.name)
        else:
            self.__dict__[attr] = val

    def __repr__(self):
        return "<%s, %s>" % (self.name,self.getinfo())

    def __getitem__(self,i):
        return self.value[i]

    def __setitem__(self,i,val):
        self.value[i] = val
        return val

    def getinfo(self,extended=False):
        "return informational string for symbol, used by __repr__"
        vtype = self.type
        if self.constant: vtype = '%s (constant)' % vtype
        sout = "%s, %s" % (vtype,repr(self.value))

        if vtype == symTypes.variable:
            nelem = datalen(self.value)
            if nelem == 1:
                sout = "%s(scalar), value=%s" % (vtype,repr(self.value))
            elif type(self.value) == Num.ArrayType:
                sout = "%s(array), npts=%i, shape=%s" % (vtype,self.value.size,self.value.shape)
            elif isinstance(self.value,(str,unicode)):
                t = str(type(self.value))[1:-1].replace('type ','').replace("'","")
                sout = "%s(%s), value=%s" % (vtype,t,repr(self.value))
            elif isinstance(self.value,(dict,tuple,list)):
                t = str(type(self.value))[1:-1].replace('type ','').replace("'","")
                sout = "%s(%s), len=%i" % (vtype,t,nelem)
        elif vtype == symTypes.group:
            sout =  "%s %i symbols, %i subgroups" % (vtype,
                                                     len(self.value.symbols()),
                                                     len(self.value.subgroups()))
        elif vtype == symTypes.pyfunc:
            sout =  "%s" % (vtype)
        elif vtype == symTypes.defvar:
            sout =  "%s, ='%s', cached value=%s" % (vtype, self.desc, repr(self.value))
            extended = False
        elif vtype == symTypes.defpro:
            args = ','.join(self.args)
            for k,v in  self.kws.items(): args = "%s,%s=%s" % (args,k,str(v))
            sout =  "%s args='%s'" % (vtype,args)
        if extended: sout = "%s %s" % (sout,self.desc)
        return sout


class SymbolTable(Group):
    def __init__(self,libs=None, writer=None, tdl=None):
        Group.__init__(self)
        self.tdl          = tdl
        self.writer       = writer  or sys.stdout

        self.loaded_libs  = []
        for i in initGroups:
            self[i] = Group(i,status='nodelete')
            
        self.initialize(libs=libs)
        self.ranlen  = 4
        
    def initialize(self,libs=None):
        self.initialized  = False

        self['_sys']['searchGroups'] = Symbol('searchGroups',value=initGroups)
        self['_sys']['LocalGroup']   = Symbol('LocalGroup', value='')
        self['_sys']['ModuleGroup']  = Symbol('ModuleGroup',value='')
        self['_sys']['path']         = Symbol('path',value=['.'])

        self.initialized = True
        
        if isinstance(libs,(tuple,list)):
            for i in libs: self.import_lib(i)
                
        
    def _getSysValue(self,name='searchGroups'):
        if len(name)<0: return None
        try:
            sys = self['_sys']
            if sys.has_key(name):
                if isSymbol(sys[name]):
                    return sys[name].value
                return sys.name
        except:
            return None

    def _getSysGroup(self,name='LocalGroup'):
        g = self._getSysValue(name)
        if self.has_key(g):
            return g
        elif self.initialized:
            self['_sys'][name].value = ''
        return None
        
    def __repr__(self):
        return "<SymbolTable id=%s>" % (hex(id(self)))

    def _lookup(self,name, isAbsolute=False, insert=None):
        """Locate and return object in symbol table given its name

        For searching, there are 'absolute' and 'find' options:

        isAbsolute=True   (not the default):
            the name must be absolute, that is referenced from the top
            the list of paths searched = [self].

        isAbsolute=False   (the default):            
            the list of paths searched =
            [self.LocalGroup,self.ModuleGroup] + self.SearchGroups + [self]

        if object is not found, there are two options

        insert = None (default):
            return None
        
        insert = 'Symbol' or 'Group' (not the default):
            create Symbol or Group object under self.LocalGroup (and all needed
            subgroups, if applicable), return object.
        """
        parts = splitName(name)
        top   = parts.pop(0)
        # print '_lookup ', name, top, parts , isAbsolute, insert

        search = []
        LocalGroup  = self._getSysGroup('LocalGroup')
        ModuleGroup = self._getSysGroup('ModuleGroup')
        # print '_lookup local/module ', LocalGroup, ModuleGroup
        
        if not isAbsolute:
            try:
                sgroups = self. _getSysValue('searchGroups')
            except:
                sgroups = initGroups
            if sgroups is None: sgroups = []
            for gname in [LocalGroup,ModuleGroup] + sgroups:
                if gname not in (None,''):
                    g = self.__get_group(name=gname)
                    if g is not None and g not in search: search.append(g)

        
        if self not in search:   search.append(self)
        # print '_lookup search groups = ', name,   search
        
        # 
        obj = None
        for s in search:
            if s.has_key(top):
                obj = s[top]
                break

        # print 'Insert? ', obj, insert
        if obj is None and insert is not None:
            obj = self
            if LocalGroup not in (None,''):
                obj = self.__get_group(name=LocalGroup)

            if len(parts) == 0:
                obj = self._createEntry(obj,top,insert)
            else:
                obj = self.__add_subgroups(obj,[top],name)

        # now, if there are sub-parts, look for these
        if obj is not None:
            for i,p in enumerate(parts):
                if obj.has_key(p):
                    obj = obj[p]
                elif insert is not None:
                    if i == len(parts)-1:
                        obj = self._createEntry(obj,p,insert)
                    else:
                        obj = self.__add_subgroups(obj,p,name)
        return obj

    def _createEntry(self,obj,name,type='Symbol'):
        if type=='Group':
            obj[name] = Group(name=name,status='normal')
        else:
            obj[name] = Symbol(name,value=None)
        return obj[name]

    def __get_group(self,obj=None,name=''):
        """ returns a group object for a (possibly nested) group name
        This requires that for a chain of attribute names (a.b.c.d.e) each
        attribute name is a valid Group

        returns None if group cannot be found
        """
        if obj is None: obj = self
        if name == '': return obj

        for p in splitName(name):
            if not obj.has_key(p): return None
            if not isGroup(obj[p]): 
                raise SymbolError, "'%s' is not a subgroup of '%s'" % (p, name)
            obj = obj[p]
        return obj
       
    def __add_subgroups(self, obj, parts, name, status='normal'):
        """ add a nested set of sub-groups for adding symbol or group

        This ensures that for a chain of attribute names (a.b.c.d.e) each
        of the names is currrently a valid Group or it creates a Group
        
        Warning: do not use for name resolution, as this will either clobber
        or miss attributes of valid python objects.
        """
        if isinstance(parts,str): parts = [parts]        
            
        for p in parts:
            if obj.has_key(p):
                if not isGroup(obj[p]): 
                    raise SymbolError, "'%s' is not a subgroup of '%s'" % (p, name)
            else:
                obj[p] = Group(name=p,status=status)
            obj = obj[p]
        return obj

    def addGroup(self,name, status='normal', isAbsolute=False):
        """ add a group:
        trying to add a group that exists will not destroy the group """
        o = self._lookup(name, isAbsolute=isAbsolute, insert='Group')
        o.status = status
        return o
    
    def placeGroup(self, group, name, status='normal'):
        """place an existing Group in the SymbolTable"""
        if not isGroup(group):
            raise SymbolError, ' cannot assign %s as a group.'% group

        parts = splitName(name)
        lastname = parts.pop()
        obj = self
        if len(parts) > 0:
            parent = '.'.join(parts)
            obj = self._lookup(parent,insert='Group')
        obj[lastname] = group
        obj[lastname].name = lastname
        
    def addTempGroup(self,prefix=None,**kw):
        " add a randomly named group, as for procedure namespaces"
        if prefix is None: prefix = '@'
        prefix.replace('.','_') # make sure it is a toplevel group!!

        ntry,nlen = 0,self.ranlen
        gname= "%s%s" % (prefix,randomName(n = nlen))
        while self.has_key(gname):    # avoid name collision!!
            gname= "%s%s" % (prefix,randomName(n = nlen))
            ntry += 1
            if (ntry > 32**nlen):         # if we have this many collisions, we're
                self.ranlen = nlen = nlen+1   # close to full, so get a larger address space
                ntry  = 0
                if (nlen > 12):
                    raise SymbolError, ' cannot create temporary group??'

        return self.addGroup(gname,status='delete')

    def clearTempGroups(self):
        "clear all toplevel groups with delete status"
        for g in self._groups():
            if g.status == 'delete': delattr(self, g)

    def getSymbol(self,name,**kw):
        return self._lookup(name,**kw)

    def getSymbolValue(self,name,**kw):
        obj = self._lookup(name,**kw)
        if isinstance(obj,Symbol):  return obj.value
        return obj
       
    def hasSymbol(self,name):
        "returns whether a Symbol exists"
        return self._lookup(name) is not None

 
    def createSymbolInGroup(self,group,name,value=None,**kw):
        """create a Symbol object inside of a Group object)"""
        if isGroup(group):
            if not group.has_key(name):
                group[name] = Symbol(name=name,value=value,**kw)
            elif value is not None:
                group[name].value = value
            return group[name]
        elif isinstance(group,(str,unicode)):
            return self.addSymbol(name = "%s.%s" % (group,name),value=value,**kw)
            

    def addSymbol(self,name, value=None,**kw):
        " add as ymbol, including passing args to created Symbol() "
        parts = splitName(name)
        suffix = parts.pop()
        if len(parts)>0:
            obj = self._lookup(".".join(parts),insert='Group')
        obj = self._lookup(name)
        
        if obj.has_key(suffix):
            print 'Warning: overwriting attribute ', suffix
        obj[suffix] = Symbol(suffix,value=value,**kw)
        return obj[suffix]

    def delSymbol(self,name):
        " delete a symbol"
        parts = splitName(name)
        suffix = parts.pop()
        if len(parts)>0:
            obj = self._lookup(".".join(parts),insert='Group')
        obj = self._lookup(name)
        if obj.has_key(suffix): obj.pop(suffix)
        return None

    def setSymbol(self,name,value,**kw):
        # print 'Symbol.setSymbol ', name, value, kw
        obj = self._lookup(name,insert='Symbol')
        # print 'set Sym lookup -> ', name, obj, isSymbol(obj)
        if isinstance(value,tuple): value = list(value)
        if isinstance(value,list):  value = list2array(value)
        
        if obj is None:
            obj = Symbol(name,value = value, **kw)
        elif isSymbol(obj):
            obj.value = value
        else:
            raise SymbolError, ' cannot set Symbol %s '% name
        
        return obj

    def setVariable(self, name, value=None, **kws):
        "add a regular variable" 
        # print 'Symbol.setVariable ', name, value, kws
        return self.setSymbol(name,value=value,type=symTypes.variable,**kws)

    def setDefVariable(self, name, expr, **kws):
        "add defined variable" 
        
        print 'Set Def Var: would compile ', expr
        val   = 'CODE: %s' % expr
        if val is None: return None
        return self.setSymbol(name,value=val, type=symTypes.defvar,code=expr,**kws)
    
    def setProcedure(self,name, code, group=None, **kws):
        if group is None:
            group = self.ModuleGroup
        return self.setSymbol(name,value=name,type=symTypes.defpro,
                              group=group, code=code, **kws)

    def setFunction(self,name,func,code=None, desc=None,cmd_out=None):
        "add a function"
        if func is None: return None
        fcn_kws = {}
        try:
            if desc is None: desc = func.__doc__
            try:
                if (func.__name__.startswith('tdl') or 
                    'tdl' in inspect.getargspec(func)[0]):
                    fcn_kws = {'kws':{'tdl':self.tdl} }
            except TypeError: # many python methods implemented in C (ie, numpy ufuncs)
                              # will raise a TypeError here.
                pass
        except:
            raise SymbolError, 'cannot add function %s ' % name

        return self.setSymbol(name,value=func,type=symTypes.pyfunc,
                              code=code,desc=desc,cmd_out=cmd_out,**fcn_kws)

    def import_lib(self,lib = None):
        " import or reload module given module name or object"
        
        if lib is None: return None

        mod, msg = None, None
        if type(lib) == types.ModuleType:
            mod = lib
            lib = lib.__name__
            
        if type(lib) == types.StringType:
            msg = '    Error loading module %s:' % lib
            try: 
                mod = __import__(lib)
                components = lib.split('.')
                for comp in components[1:]:
                    mod = getattr(mod, comp)
            except ImportError:
                msg = '    Error loading module %s:' % lib
        if mod in self.loaded_libs:
            try:
                mod = reload(mod)
            except ImportError:
                msg = '    Error loading module %s:' % lib
        if mod is None:
            self.writer.write("    cannot load module %s !" % lib)
            if msg is not None: self.tdl.ShowError(msg, showtraceback=False)
            return None

        # mod is now a real module, not a string of the module name
        self.writer.write("    loading %s ..." % mod.__name__)
        self.writer.flush()
        if mod not in self.loaded_libs: self.loaded_libs.append(mod)
        
        try:
            if hasattr(mod,'_groups_'):
                for nam,toplevel in getattr(mod,'_groups_',(None,True)):
                    if nam: self.addGroup(nam,toplevel=toplevel)
            if hasattr(mod,'_var_'):
                for nam,val in getattr(mod,'_var_',{}).items():
                    self.setSymbol(nam,val)
            if hasattr(mod,'_consts_'):
                for nam,val in getattr(mod,'_consts_',{}).items():
                    self.setSymbol(nam,val,constant=True)
            if hasattr(mod,'_func_'):
                for nam,val in getattr(mod,'_func_',{}).items():
                    cmdOut = None
                    func   = val
                    if type(val) == types.TupleType:
                        func = val[0]
                        if len(val) > 1: cmdOut = val[1]
                    x =self.setFunction(nam,func,cmd_out=cmdOut) 
            if hasattr(mod,'_scripts_'):
                for nam in getattr(mod,'_scripts_',[]):
                    try:
                        file_path = os.path.abspath(os.path.dirname(mod.__file__))
                        file_name = os.path.join(file_path,nam)
                        if os.path.exists(file_name) and os.path.isfile(file_name):
                            self.tdl.load_file(file_name)
                        else:
                            self.writer("Warning: Cannot find lib script file: %s\n" % file_name)
                    except:
                        self.tdl.ShowError(msg="Error loading script file '%s'"  % file_name)
            if hasattr(mod,'_help_'):
                for nam,val in getattr(mod,'_help_',{}).items():
                    self.tdl.help.add_topic(nam,val)
            if hasattr(mod,'_init_'):
                f = getattr(mod,'_init_')
                fname = '_@@%s@@init_%s' % (mod.__name__, randomName())
                fname = fname.replace('.','_')
                s = self.setFunction(fname, f)
                self.tdl.eval("%s()" % fname)
                self.delSymbol(fname)
            import_msg = '   ok.'
        except ImportError:
            import_msg = '   import failed!'

        self.writer.write(" %s\n" % import_msg)


#     def import_lib(self,lib):
#         " import or reload module given module name or object"
# 
#         if lib is None: return None
#         if self.tdl is None: return None
# 
#         mod, msg = None, None
#         if type(lib) == types.ModuleType:
#             mod = lib
#             lib = lib.__name__
#             
#         if type(lib) == types.StringType:
#             msg = '    Error loading module %s:' % lib
#             try: 
#                 mod = __import__(lib)
#                 components = lib.split('.')
#                 for comp in components[1:]:
#                     mod = getattr(mod, comp)
#             except ImportError:
#                 msg = '    Error loading module %s:' % lib
#         if mod in self.loaded_libs:
#             try:
#                 mod = reload(mod)
#             except ImportError:
#                 msg = '    Error loading module %s:' % lib
#         if mod is None:
#             self.writer.write("    cannot load module %s !" % lib)
#             if msg is not None: self.tdl.ShowError(msg, showtraceback=False)
#             return None
# 
#         # mod is now a real module, not a string of the module name
#         title = getattr(mod,'title',mod.__name__)
#         self.writer.write("    loading %s ..." % title)
#         self.writer.flush()
#         if mod not in self.loaded_libs: self.loaded_libs.append(mod)
#         
#         try:
#             if hasattr(mod,'_groups_'):
#                 for nam,toplevel in getattr(mod,'_groups_',(None,True)):
#                     if nam: self.addGroup(nam,toplevel=toplevel)
#             if hasattr(mod,'_var_'):
#                 for nam,val in getattr(mod,'_var_',{}).items():
#                     self.setSymbol(nam,val)
#             if hasattr(mod,'_consts_'):
#                 for nam,val in getattr(mod,'_consts_',{}).items():
#                     self.setSymbol(nam,val,constant=True)
#             if hasattr(mod,'_func_'):
#                 for nam,val in getattr(mod,'_func_',{}).items():
#                     cmdOut = None
#                     func   = val
#                     if type(val) == types.TupleType:
#                         func = val[0]
#                         if len(val) > 1: cmdOut = val[1]
# 
#                     x =self.setFunction(nam,func,cmd_out=cmdOut) 
#             if hasattr(mod,'_scripts_'):
#                 for nam in getattr(mod,'_scripts_',[]):
#                     try:
#                         file_path = os.path.abspath(os.path.dirname(mod.__file__))
#                         file_name = os.path.join(file_path,nam)
#                         if os.path.exists(file_name) and os.path.isfile(file_name):
#                             self.tdl.load_file(file_name)
#                         else:
#                             self.writer("Warning: Cannot find lib script file: %s\n" % file_name)
#                     except:
#                         self.tdl.ShowError(msg="Error loading script file '%s'"  % file_name)
#             if hasattr(mod,'_help_'):
#                 for nam,val in getattr(mod,'_help_',{}).items():
#                     self.tdl.help.add_topic(nam,val)
#             if hasattr(mod,'_init_'):
#                 f = getattr(mod,'_init_')
#                 fname = '_@@%s@@init_%s' % (mod.__name__, randomName())
#                 fname = fname.replace('.','_')
#                 s = self.setFunction(fname, f)
#                 self.tdl.eval("%s()" % fname)
#                 self.delSymbol(fname)
#             import_msg = 'ok.'
#         except ImportError:
#             import_msg = 'import failed!'
# 
#         self.writer.write(" %s\n" % import_msg)
#         


#     def initialize(self,libs=None):
#         ll = ('TdlBuiltins','TdlNumLib')
#         if libs is not None and type(libs) in (types.ListType,types.TupleType):
#             ll = ll + list(libs)
#         # 
#         spath = self.getSymbolValue('_sys.path')
# 
# 
#         if type(spath) == types.ListType:
#             for p in spath: set_path(p)
#             
#         elif type(spath) == types.StringType:
#             set_path(spath)
#         # import all libraries
#         for i in ll:
#             self.import_lib(i)
        
        


