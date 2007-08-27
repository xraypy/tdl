# M. Newville Univ of Chicago (2006)
#
# -------------
# Modifications
# -------------
#
# Sept 10, 2006 MN  complete rewrite (version 0.3)
#
#########################################################################

__version__  = '0.6'

# symbol Groups that must always be present, and starting search order
initGroups    = ('_main','_sys','_math','_plot','_builtin')
init_loadLibs = ['TdlBuiltins','TdlNumLib']
_TopGroupName = '_TDL_'

import types
import os
import random
import sys
import inspect
import re
from copy import deepcopy

import Num
from Num import num_version, ndarray
from Util import find_unquoted_char, split_delim, datalen, set_path, list2array
from Util import PrintExceptErr, SymbolError, ConstantError, show_more
from string import ascii_lowercase, digits

isValidName = re.compile(r'^[a-zA-Z_\$&@][a-zA-Z_\$&@0-9]*$').match
def randomName(n=6):
    return ''.join([(ascii_lowercase + digits)[random.randrange(36)] for i in range(n)])

class symTypes:
    variable  = 'variable'
    group     = 'symbol group'
    pyobj     = 'python object'
    file      = 'file'
    defvar    = 'defined variable'
    pyfunc    = 'python function'
    defpro    = 'tdl procedure'
    Data      = (variable,group,pyobj,file,defvar)
    Funcs     = (pyfunc,defpro)
    All       = Data + Funcs
    NonGroups = (variable,pyobj,file,defvar,pyfunc,defpro)


var_typecodes = {types.StringType:'string',
                 types.IntType: 'int',
                 types.LongType:'int',
                 types.FloatType:'float',
                 types.ComplexType:'complex',
                 types.ListType:'list',
                 types.DictType:'dict',
                 ndarray:'array'}

class Group(dict):
    """ symbol group is an extended dictionary that holds symbols and other symbol groups"""
    type = symTypes.group
    def __init__(self,name=None,filename=None,toplevel=False,status='normal', vars=None):
        self.filename = filename
        self.status   = status
        self.toplevel = toplevel
        self.setname(name)
        if vars is not None:
            for k,v in vars.items():
                # Remove automatic typecast of lists as arrays
                #if type(v) == types.ListType: v = list2array(v)
                self.setSymbol(k,v)
        
    def setname(self,name=None):
        if name is None: name = ''
        self.name = name

    def __repr__(self):
        return "<%s, %s>" % (self.name, self.getinfo())
    
    def getinfo(self,extended=False):
        # return "status=%s, id=%s" % (self.status, hex(id(self)))
        s = "group(%s)" % self.status
        if extended:
            nv,nf,ng = self.stats()            
            s = "%s: %i variables, %i functions, %i subgroups" % (s,nv,nf,ng)
        return s

    def addGroup(self,name,group=None,filename=None,status='normal',toplevel=False):
        if self.status=='frozen' and not self.has_key(name):
            raise SymbolError, ' group %s is frozen, and cannot be extended.' % self.name
        fn = filename or self.filename
        sname = "%s.%s" %(self.name,name)
        if toplevel: sname = name
        if isGroup(group):
            self[name] = group
            self[name].setname(name)
        else:
            self[name] = Group(name=sname,filename=fn,status=status)
        return self[name]

    def delGroup(self,name):
        # print 'delGroup ', self.name, name, self.has_key(name), type(self[name])
        if self.status=='frozen': raise SymbolError, ' group %s is frozen.' % self.name
        if (self.has_key(name) and isGroup(self[name]) and self[name].status != 'nodelete'):
            self.pop(name)
        else:
            raise SymbolError, ' cannot delete group %s.%s' % (self.name, name)

    def setSymbol(self,name,value=None,**kw):
        # print 'Group: setSymbol ', self.name, name, value, kw
        if self.status=='frozen': raise SymbolError, ' group %s is frozen.' % self.name
        if (self.has_key(name) and isSymbol(self[name]) and self[name].constant):
            raise ConstantError, ' cannot overwrite constant %s' % name
        self[name] = Symbol(name="%s.%s" %(self.name,name),value=value,**kw)
        return self[name]
        
    def Symbols(self): return self.keys()
    
    def hasSymbol(self,name): return self.has_key(name)

    def delSymbol(self,name):
        if self.status=='frozen': raise SymbolError, ' group %s is frozen.' % self.name
        if self.has_key(name) and isSymbol(self[name]):
            self.pop(name)

    def stats(self):
        "return (n_variables, n_functions, n_groups) in a group"
        nvar, nfunc, ngroup = 0,0,0
        for sym in self.values():
            if  isGroup(sym):
                ngroup = ngroup + 1
            elif sym.type in (symTypes.defpro,symTypes.pyfunc):
                nfunc  = nfunc  + 1
            else: nvar = nvar + 1
        return (nvar,nfunc, ngroup)


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
       mod        module namespace for defpro
       cmd_out    method for processing output when a pyfunc of 
                     defpro is run as a command.

    Notes:
       a Simple Name matches [a-zA-Z_\$][a-zA-Z_\$0-9]*
       Of course, a Symbol may be part of a hiearchy of Groups,
       and so may have a much longer Full Name.
    """
    name     = ''
    constant = False
    def __init__(self, name, value=None, constant= False,
                 stype=symTypes.variable,
                 desc=None, help=None, cmd_out=None,
                 code=None, args=None, mod=None, kws=None):
        self.constant = False
        self.name     = name
        self.cmd_out  = cmd_out
        self.value    = value
        self.code     = code
        self.desc     = desc or ''
        self.help     = help or ''
        self.args     = args or []
        self.kws      = kws  or {}
        self.mod      = mod  or '_main'
        
        # check that type is valid
        if stype not in symTypes.All:
            raise SymbolError, 'cannot add symbol "%s" with type=%s' % (name,stype)
        self.type = stype
        self.constant = constant

    def getHelp(self):
        return self.help

    def getCode(self):
        return deepcopy(self.code)

    def call(self, *args,**kws):
        if self.type != symTypes.pyfunc: return self.value

        x = self.kws.copy() ; x.update(kws)
        val = self.value(*args,**x)
        return val

    def __call__(self, *args,**kws): self.call(*args,**kws)
    
    def cmdout(self,val,**kws):
        if val is None: return None
        if self.cmd_out is None: return str(val)

        x = self.kws.copy() ; x.update(kws)
        return self.cmd_out(val,**x)
         
    def __setattr__(self, attr, val):
        """ here to prevent re-setting of constants"""
        if self.constant and attr == 'value':
            raise ConstantError,  'cannot set value of constant %s' % (self.name)
        self.__dict__[attr] = val

    def __repr__(self):
        return "<%s, %s>" % (self.name,self.getinfo())
    
    def getinfo(self,extended=False):
        "return informational string for symbol, used by __repr__"
        
        vtype = self.type
        if self.constant: vtype = '%s (constant)' % vtype
        #sout = "%s, %s" % (vtype,repr(self.value))
        sout = "%s" % (vtype)

        if vtype == symTypes.variable:
            ty = type(self.value)
            if (ty in var_typecodes.keys()) or (ty in Num.typeDict.values()): 
                nelem = datalen(self.value)
                if (nelem == 1): 
                    sout = "%s(scalar), value=%s" % (vtype,repr(self.value))
                elif type(self.value) == Num.ArrayType:
                    sout = "%s(array), npts=%i, shape=%s" % (vtype,self.value.size,self.value.shape)
                elif type(self.value) == types.StringType:
                    sout = "%s(string), value=%s" % (vtype,repr(self.value))
                elif type(self.value) in (types.DictType,types.TupleType,types.ListType):
                    t = str(type(self.value))[1:-1].replace('type','').strip()
                    t = t.replace("'","")
                    sout = "%s(%s), len=%i" % (vtype,t,nelem)
            else:
                sout = "%s (%s)" % (vtype, ty)
        elif vtype == symTypes.pyfunc:
            sout =  "%s" % (vtype)
        elif vtype == symTypes.defvar:
            sout =  "%s, ='%s', cached value=%s" % (vtype, self.desc, repr(self.value))
            extended = False
        elif vtype == symTypes.defpro:
            args = ','.join(self.args)
            for k,v in  self.kws.items(): args = "%s,%s=%s" % (args,k,str(v))
            sout =  "%s args='%s'" % (vtype,args)

        if extended:
            sout = "%s\n  %s" % (sout,self.desc)

        return sout

############################################################

def splitname(s):
    parts = s.split('.')
    for i,p in enumerate(parts):
        if i == 0 and p== _TopGroupName:
            pass
        elif not isValidName(p):
            msg = 'invalid name "%s"' % (s)
            if len(s)>1: msg = 'invalid name "%s" (error at "%s")' % (s,p)
            raise SymbolError, msg
    return parts

def _splitName(name,check=True):
    parts = name.split('.')
    if check:
        for i in parts:
            if not isValidName(i):
                msg = 'invalid name "%s"' % (name)
                if len(n)>1: msg = 'invalid name "%s" (error at "%s")' % (name,i)
                raise SymbolError, msg
    return parts

def isGroup(x):
    try:
        return (x.type == symTypes.group)
    except:
        return False
    
def isSymbol(x):
    return (hasattr(x,'value') and hasattr(x,'getCode'))

class SymbolTable:
    def __init__(self,libs=None, writer=None, tdl=None, init=True):
        self.tdl          = tdl
        self.writer       = writer  or sys.stdout
        self.loaded_libs  = []
        self.searchGroups = []

        self.data = Group(name=_TopGroupName, status='nodelete')

        self.LocalGroup  = '_main' 
        self.ModuleGroup = '_main' 
        for i in initGroups:
            self.addGroup(i,toplevel=True)
            self.searchGroups.append(i)

        self.setSymbol('_sys.path',['.'])
        self.setSymbol('_sys.searchGroups',self.searchGroups)
        
        if init: self.initialize()

    def initialize(self,libs=None):
        ll = ('TdlBuiltins','TdlNumLib')
        if libs is not None and type(libs) in (types.ListType,types.TupleType):
            ll = ll + list(libs)
        # 
        spath = self.getSymbolValue('_sys.path')

        if type(spath) == types.ListType:
            for p in spath: set_path(p)
        elif type(spath) == types.StringType:
            set_path(spath)
        # import all libraries
        for i in ll:
            self.import_lib(i)

    def import_lib(self,lib):
        " import or reload module given module name or object"

        if self.tdl is None: return None

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
        title = getattr(mod,'title',mod.__name__)
        self.writer.write("    loading %s ..." % title)
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

    def reimport_libs(self):
        """Reload all imported libraries"""
        for mod in self.loaded_libs:
            name = mod.__repr__()
            name = name[1:name.find('from')]
            try:
                print "Reimport: %s" % name
                self.import_lib(mod)
            except:
                print "\n    **Failed to reimport %s\n" % name
        return

    def _lookup(self,name):
        """ simple lookup of where a name in the symbol table is:
        """
        parts = splitname(name)
        p     = parts[0]
        parent,sym = None,None

        if (self.data.has_key(self.LocalGroup) and 
            self.data[self.LocalGroup].has_key(p)):
            parent = self.data[self.LocalGroup]
        elif (self.LocalGroup != self.ModuleGroup and
              self.data.has_key(self.ModuleGroup) and 
              self.data[self.ModuleGroup].has_key(p)):
            parent = self.data[self.ModuleGroup]            
        elif p == _TopGroupName:
            parent,sym = None,self.data
        elif p in self.data.keys():
            parent = self.data
        else:
            for g in self.searchGroups:
                if self.data[g].has_key(p):
                    parent = self.data[g]
                    break

        if parent is not None: sym = parent[p]
        if sym is not None: parts.pop(0)

        while len(parts)>0:
            p = parts[0]
            if (isGroup(sym) and sym.has_key(p) and
                 (isGroup(sym[p]) or isSymbol(sym[p]))):
                parent,sym= sym,sym[p]
            else:
                break
            parts.pop(0)
        return (parent,sym,parts)

    def _normalize_sym(self,sym,toplevel=False):
        " a common task: look up a group name or select the LocalGroup"
        if sym is None:
            if toplevel or not self.data.has_key(self.LocalGroup):
                sym = self.data
            else:
                sym = self.data[self.LocalGroup]
        return sym
        
    def addGroup(self, name, toplevel=False, status='normal'):
        parent,sym,parts = self._lookup(name)
        sym = self._normalize_sym(sym,toplevel=toplevel)
        if not isGroup(sym):
            raise SymbolError, ' cannot create group %s '% name
        for p in parts:  sym = sym.addGroup(p,status=status,toplevel=toplevel)
        return sym

    def clearTempGroups(self):
        " clear all toplevel groups with delete status"
        k = self.data.keys()
        for i in k:
            if type(self.data[i]) == Group and self.data[i].status=='delete':
                self.data.pop(i)
    
    def addTempGroup(self,prefix=None,nlen=6,**kw):
        " add a randomly named group, as for procedure namespaces"
        if prefix is None: prefix = ''
        if nlen < 3: nlen  = 3
        gname= "@%s_%s" % (prefix,randomName(n = nlen))
        ntry = 0
        while self.data.has_key(gname):  # avoid name collision!!
            gname= "@%s_%s" % (prefix,randomName(n = nlen))
            ntry = ntry+1
            if (ntry > 10**nlen):  nlen = nlen + 1
            
        return self.addGroup(gname,toplevel=True,status='delete')

    def placeGroup(self, group, name, toplevel=False,status='normal'):
        """ place an existing Group in the SymbolTable"""
        if not isGroup(group):  raise SymbolError, ' cannot assign %s as a group.'% group

        parent,sym,parts = self._lookup(name)

        nam = sym.name.split('.')[-1]
        parent.delSymbol(name)
        parent.addGroup(nam,group=group,status=status)
        parent[nam].setname(nam)
        for elem in  parent[nam].values():
            names = elem.name.split('.')
            elem.name = "%s.%s" % (nam,names[-1])

    def __gather(self,x, _type=None,recurse=True):
        for k,v in x.items():
            if isSymbol(v):
                if v.type in _type: self.workbuffer.append(k)
            elif isGroup(v):
                self.workbuffer.append(v.name)
                if recurse: self.__gather(v,_type=_type)

    def __list(self,_type=symTypes.All):
        self.workbuffer = []
        self.__gather(self.data,_type=_type)
        return self.workbuffer

    def listAll(self):        return self.__list(_type=symTypes.All)
    def listGroups(self):     return self.__list(_type=symTypes.group)
    def listVariables(self):  return self.__list(_type=symTypes.Data)
    def listFunctions(self):  return self.__list(_type=symTypes.Funcs)
        
    def getStats(self,groupname=None):
        'returns statistics about toplevel symbol table'
        out = {}

        def groupcount(g):
            nvar, nfunc, ngroup = 0,0,0
            for sym in g.values():
                if  isGroup(sym):
                    ngroup = ngroup + 1
                elif sym.type in (symTypes.defpro,symTypes.pyfunc):
                    nfunc  = nfunc  + 1
                else: nvar = nvar + 1
            return (nvar,nfunc, ngroup)
            
        for sym in self.data.values():
            if isGroup(sym):
                nam = sym.name
                if nam.startswith('%s.' % _TopGroupName):
                    nam = nam[len(_TopGroupName)+1:]
                out[nam] = groupcount(sym)
            #elif sym.type in (symTypes.defpro,symTypes.pyfunc):
            #    nfunc  = nfunc  + 1
            #else: nvar = nvar + 1
        return out
            
    def setSymbol(self,name,value,create=True,**kw):
        parent,sym,parts = self._lookup(name)
        sym = self._normalize_sym(sym)
        
        if isSymbol(sym):
            sym.value = value
            for k,v in kw.items():
                if hasattr(sym,k): setattr(sym,k,v)

        elif isGroup(sym) and len(parts)>0:
            if len(parts) > 1:
                if create:
                    while len(parts)>1:
                        p = parts.pop(0)
                        sym = sym.addGroup(p)
                else:
                    raise SymbolError, 'cannot add Symbol %s: would need to create subgroups.' % (name,sym.name)
            sym = sym.setSymbol(parts[0],value,**kw)
        else:
            raise SymbolError, ' cannot add Symbol  %s: %s is not a group '% (name, sym.name)
        return sym
        
    def setVariable(self, name, value=None, **kws):
        "add a regular variable" 
        return self.setSymbol(name,value=value,stype=symTypes.variable,**kws)
    
    def setDefVariable(self, name, code, desc,**kws):
        "add defined variable" 
        try:
            val = self.tdl.expr_eval(code)
        except:
            val = None        
        return self.setSymbol(name,value=val, stype=symTypes.defvar,
                              code=code,desc=desc,**kws)

    def setProcedure(self,name,code,desc=None,**kws):
        return self.setSymbol(name,value=name,stype=symTypes.defpro,
                              mod=self.ModuleGroup,
                              code=code, desc=desc,**kws)

    def setFunction(self,name,func,ftype=symTypes.pyfunc,code=None,
                    desc=None,cmd_out=None):
        "add a function"
        if func is None: return None

        if ftype not in symTypes.Funcs:
            raise SymbolError, 'cannot add function %s with type  %s' % (name,ftype)

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
        return self.setSymbol(name,value=func,stype=ftype,code=code,
                              desc=desc,cmd_out=cmd_out,**fcn_kws)

    def getSymbolLocalGroup(self,name):
        parts = splitname(name)
        p     = parts[0]
        parent,sym = None,None

        if self.data.has_key(p):
            sym = self.data[p]
            parts.pop(0)
        elif self.data.has_key(self.LocalGroup):
            sym = self.data[self.LocalGroup]
        else:
            sym = self.data
            
        if isGroup(sym) and len(parts)>1:
            while len(parts)>1:
                p = parts.pop(0)
                if sym.has_key(p):
                    sym = sym[p]
                else:
                    sym = sym.addGroup(p)

        if sym.hasSymbol(parts[0]):return sym[parts[0]]

        return sym.setSymbol(parts[0],None)

    def getSymbol(self,name):
        parent,sym, parts = self._lookup(name)
        # print 'getSymbol ', name, self.LocalGroup,self.ModuleGroup
        # print '   --> ', parent, sym, parts
        if sym is None:   sym = self._normalize_sym(sym)

        if len(parts)==1 and isGroup(sym):
            if sym.has_key(parts[0]):  return sym[parts[0]]

        if sym is  None: sym  = self.data[self.LocalGroup]

        # print '   Sym  / Parts ', sym, parts        
        for i in parts:
            stmp  = None
            try:
                stmp = getattr(sym.value,i)
            except:
                try:
                    stmp = getattr(sym,i)
                except:
                    pass
            if stmp is None:
                raise SymbolError, 'cannot find or create "%s" in "%s"' % (i,sym)
            sym = stmp
        return sym
            
    def getSymbolValue(self,name):
        n = self.getSymbol(name)
        if n is None:   return None
        if isSymbol(n): return n.value
        return n

    def getVariable(self, name):
        return self.getSymbol(name)

    def getFunction(self,name):
        "returns a function, assuming it exists"
        ret = None
        try:
            f = self.getSymbol(name)
            if isSymbol(f):
                if f.type in symTypes.Funcs: ret = f 
            elif callable(f):
                ret = f
        except:
            pass
        return ret

    def getGroup(self,name):
        parent,sym,parts = self._lookup(name)
        if not isGroup(sym) or parts != []:
            raise SymbolError, ' cannot get group %s '% name
        return sym

    def hasSymbol(self,name):
        "returns whether a Symbol exists"
        try:
            return self.getSymbol(name) is not None
        except SymbolError:
            return False

    def hasVariable(self,name):
        return self.hasSymbol(name)
    
    def hasFunction(self,name):
        "returns whether a function exists"
        try:
            return self.getFunction(name) is not None
        except SymbolError:
            return False

    def hasGroup(self,name):
        "returns whether a function exists"
        try:
            return self.getGroup(name) is not None
        except SymbolError:
            return False

    def delSymbol(self,name):
        "delete a symbol"
        parent,sym,parts = self._lookup(name)
        name = splitname(name).pop()

        # print 'del Symbol ', name, sym, parent, isGroup(sym), isSymbol(sym)
        if not isGroup(parent):
            raise SymbolError, 'cannot resolve symbol for delete %s' % name
        if isGroup(sym):
            if sym.status  == 'nodelete':
                raise SymbolError, 'cannot delete group %s  (nodelete)' % name            
            parent.delGroup(name)
        elif isSymbol(sym):
            if len(parts)>0:
                raise SymbolError, 'cannot delete symbol: "%s" not resolved (python attribute??)' % name
            parent.delSymbol(name)
            
    def delGroup(self,name):
        self.delSymbol(name)
        
    def _groupstats(self,group):
        nvar, nfunc, ngroup = 0,0,0
        for sym in group.values():
            if  isGroup(sym):
                ngroup = ngroup + 1
            elif sym.type in (symTypes.defpro,symTypes.pyfunc):
                nfunc  = nfunc  + 1
            else: nvar = nvar + 1
        return (nvar,nfunc, ngroup)
    
    def showTable(self,skip=None,do_print=True):
        print '::Symbol Table::'
        ignoregroups = []
        if skip is not None: ignoregroups = skip
        indent = 1
        buff = []
        def _showGroup(grp,indent=1, show_subgroups=False):
            st    = self._groupstats(grp)
            otab  = '   '*indent
            vtab  = '   '*(indent+1)
            nam   = grp.name + ' '*(16-len(grp.name))
            buff.append('%s%s: %i variables, %i functions, %i subgroups' % (otab,nam,st[0],st[1],st[2]))
            if show_subgroups:
                for k,v in grp.items():
                    if isGroup(v):
                        if v not in ignoregroups:  _showGroup(v,indent+1)
                    elif isSymbol(v):
                        buff.append("%s%s = %s\n" % (vtab, k, v.getinfo()))

        for g in self.data.values():  _showGroup(g,indent,show_subgroups=True)
        print '=================='
        if do_print:
            show_more(buff)
        return buff

