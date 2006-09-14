# M. Newville Univ of Chicago (2006)
#
# -------------
# Modifications
# -------------
#
# Sept 10, 2006 MN  complete rewrite (version 0.3)

#
# 6-4-06 T2
# In SymbolTable.import_lib added a case for _scripts_
# If the lib module defines a list _scripts_ = [file1,file2]
# then these will be loaded (eg to define various defs or data..)
# Its assumed that these scripts live in the same directory as the
# library module.  
#
# 4-2-06 T2
# Modified symbolTable.initialize and symbolTable.import_lib
# to correctly handle strings as module names
#
# * 2-12-06 T2
# Added new methods for getting/putting to symbol table
# back to a single dictionary of symbols
#
#########################################################################

# symbol Groups that must always be present, and starting search order
initGroups    = ('_sys','_math','_main','_plot','_builtin')

init_loadLibs = ['TdlBuiltins','TdlNumLib']


import types
import os
import random
import sys
import inspect
import re
from copy import deepcopy

from Num import Num
from Util import find_unquoted_char, split_delim, datalen, set_path
from Util import PrintExceptErr, PrintShortExcept, SymbolError, ConstantError
from string import ascii_lowercase, digits

__version__  = '0.3'

isValidName = re.compile(r'^[a-zA-Z_\$&@][a-zA-Z_\$&@0-9]*$').match

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


class symGroup(dict):
    """ variable group
    """
    type = 'group'
    def __init__(self,name=None,filename=None,toplevel=False,status='normal'):
        self.filename = filename
        self.status   = status
        self.toplevel = toplevel
        self.setname(name)

    def setname(self,name=None):
        if name is not None: self.name     = name
        self.value    = "<symGroup %s status=%s, id=%s>" % (self.name,
                                                            self.status,
                                                            hex(id(self)))
    def __repr__(self):
        self.setname(self.name)
        return self.value
        
    def addGroup(self,name,filename=None,status=None,toplevel=False):
        if self.status=='frozen' and not self.has_key(name):
            raise SymbolError, ' group %s is frozen, and cannot be extended.' % self.name
        fn = filename or self.filename
        st = status   or self.status
        sname = "%s.%s" %(self.name,name)
        if toplevel: sname = name
        self[name] = symGroup(name=sname,filename=fn,status=st)
        
    def delGroup(self,name):
        dodel = False
        if self.has_key(name) and type(self[name]) == symGroup:
            if self[name].status != 'nodelete': dodel = True
        if dodel:  self.pop(name)

    def addSymbol(self,name,value=None,**kw):
        if self.status=='frozen' and not self.has_key(name):
            raise SymbolError, ' group %s is frozen, and cannot be extended.' % self.name
        if self.has_key(name):
            if isinstance(self[name],Symbol) and self[name].constant:
                raise ConstantError, ' cannot overwrite constant %s' % name
        self[name] = Symbol(name="%s.%s" %(self.name,name),value=value,**kw)

    def delSymbol(self,name):
        if self.has_key(name) and isinstance(self[name],Symbol):
            self.pop(name)

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
                 type=symTypes.variable,
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
        if type not in symTypes.All:
            raise SymbolError, 'cannot add symbol "%s" with type=%s' % (name,type)
        self.type = type
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
        name = self.name
        if self.type == symTypes.variable:
            nelem = datalen(self.value)
            vtype = self.type
            if self.constant: vtype = 'variable (constant)'
            if nelem == 1:
                return "<%s %s: type=scalar, value=%s>" % (vtype,name,repr(self.value))
            elif type(self.value) == Num.ArrayType:
                return "<%s %s: type=array, npts=%i, shape=%s>" % (vtype,
                                                                   name,
                                                                   self.value.size,
                                                                   self.value.shape)
            else:
                t = str(type(self.value))[1:-1].replace('type','')
                t.strip()
                return "<%s %s: type=%s, len=%i>" % (vtype,name,t,nelem)
        elif self.type == symTypes.defvar:
            return "<Defined Variable %s ='%s'>" % (name,self.desc)
        elif self.type == symTypes.pyfunc:
            vtype = 'Function'
            if self.constant: vtype = 'Function (constant)'
            return "<%s %s>" % (vtype, name)
        elif self.type == symTypes.defpro:
            args = ','.join(self.args)
            for k,v in  self.kws.items(): args = "%s,%s=%s" % (args,k,str(v))
            return "<Procedure %s: args='%s'>" % (name,args)
        return "<Symbol %s: %s : %s>" % (name, self.type,repr(self.value))


############################################################
class SymbolTable:
    """
    Table of symbols and groups storing all functions and variables in tdl
    """
    vchars = ascii_lowercase + '_' + digits
    rrange = random.randrange
    def randomName(self,n=8):
        return "_%s" % ''.join([self.vchars[self.rrange(0,37)] for i in range(n)])


    def __init__(self,libs=None,writer=None,tdl=None,init=True,**kws):
        self.tdl          = tdl
        self.writer       = writer  or sys.stdout
        self.loaded_libs  = []
        self.data         = symGroup(name='_top',status='nodelete')
        self.searchGroups = []
        for i in initGroups:
            self.data.addGroup(i,status='nodelete',toplevel=True)
            self.searchGroups.append(i)
        self.LocalGroup   = '_main'
        self.ModuleGroup  = '_main'

        self.addSymbol('_sys.path',['.'])
        self.addSymbol('_sys.searchGroups',self.searchGroups)
        # usually, we're ready to initialize the SymbolTable right away
        if init:  self.initialize()
        
    def __splitName(self,name):
        n = name.split('.')
        for i in n:
           if not isValidName(i):
               raise SymbolError, 'invalid name "%s"' % (name)
        return n

    def resolveGroup(self, name):
        "returns fully resolved group object "
        names = self.__splitName(name)
        grp  = self.data.get(names[0])
        # print 'resolvegroup name ' , grp, name, names       
        if len(names)>1:
            names = names[1:]
            gn = grp.get(names[0])
            if gn is not None and type(gn)==symGroup:
                grp = gn
                names = names[1:]
            else:
                raise SymbolError, "group %s not resolved" % name
        return grp

    def __findName(self,name):
        """find symbol name for adding Symbol of Group
           Symbol Names are   ZeroOrMore(SimpleName + '.') + SimpleName

           if name is a SimpleName (no '.'): add to LocalGroup

           if name contains multiple parts (x.y.z):
               assert that the first name is in LocalGroup
               if the first name is found as a toplevel group, use that
               loop through rest of the name parts, making sure that
               all subgroups are properly found.

           returns symGroup and member name.  Examples:
                in                    out
                a                 (self.data[self.LocalGroup], 'a')
                _main.a           (self.data['_main'], 'a')
                _main.a           (self.data['_main'], 'a')
                _main.x.a         (self.data['_main']['x'], 'a')

           for the last example, if self.data['_main']['x'] does not
           exist, or is not a symGroup, a SymbolError is raised.
        
        """
        grp = self.resolveGroup(self.LocalGroup)
        names = self.__splitName(name)
        if len(names)>1 and self.data.has_key(names[0]):
            grp = self.data[names[0]]
            names = names[1:]
        while len(names) > 1:
            gn = grp.get(names[0])
            if isinstance(gn,Symbol): gn = gn.value
            if gn is None:
                raise SymbolError, 'cannot find group "%s" in "%s"' % (names[0],grp.name)
            if type(gn) != symGroup:
                print gn, type(gn)
                raise SymbolError, 'symbol "%s.%s" is not a group' % (grp.name, names[0])
            grp = gn
            names = names[1:]
        return (grp,names[0])


    def resolveName(self,name):
        """ returns  'full path to symbol', so that the symbol
        can easily be resolved from it's split absolute name.
        """
        names = self.__splitName(name)
        ret  = []
        # first look in localgroup and modulegroup - most likely case        
        if len(names)==1:
            for i in (self.LocalGroup,self.ModuleGroup):
                try:
                    sym = self.resolveGroup(i).get(name)
                    if sym is not None:
                        ret = self.__splitName(i)
                        ret.append(name)
                        return (ret,[])
                except AttributeError:
                    pass
        # now construct the search groups, and look through them
        # for the FirstName
        slist = []
        for i in ([self.LocalGroup,self.ModuleGroup] +
                  self.searchGroups + ['_sys','_builtin','_main']):
            g = self.resolveGroup(i)
            r = i,g
            if r not in slist: slist.append(r)
        # add toplevel group
        slist.append(('_top',self.data))

        for i in self.data.keys():
            g = self.resolveGroup(i)
            r = i,g
            if r not in slist: slist.append(r)

        #
        # For SimpleNames and FirstName of CompoundNames:
        #     search all searchGroups for name
        sym = None
        name0 = names[0]
        names = names[1:]
        for gname,group in slist:
            if name0 == group.name:  # name0 is the name of a group in searchGroup
                sym = group
                ret = [name0]
                break
        if sym is None:
            for gname,group in slist:
                sym = group.get(name0)   # name0 is contained in a group in searchGroup
                if sym is not None:
                    ret = [group.name,name0]
                    if group.name=='_top':  ret = [name0]
                    break
        if sym is None:
            raise SymbolError, 'could not resolve name "%s"' % (name)

        #
        # Continue with compound Names
        names.reverse()
        nx   = []
        rest = []
        for i in range(len(names)):
            n  = names.pop()
            sn = None
            if type(sym)==symGroup:
                sn = sym.get(n)
            else:
                rest = [n] + names
                break
            if sn is None:
                msg = sym.name
                raise SymbolError, 'cannot find "%s" in "%s"' % (n,msg)
            sym = sn
            nx.append(n)
        ret = ret + nx
        return (ret,rest)

    def addSymbol(self, name, value=None,**kws):
        """ basic form for adding new Symbols into Symbol Table"""
        g,n = self.__findName(name)
        g.addSymbol(name=n,value=value,**kws)
        return g[n]


    def getSymbolLocalGroup(self, name):
        """ basic form for adding new Symbols into Symbol Table"""
        g,n = self.__findName(name)
        sname = "%s.%s" % (g.name,n)
        if not self.hasSymbol(sname):
            g.addSymbol(name=n)
            s = g[n]
        else:
            s = self.getSymbol(sname)
        # print 'gslp ' , s  
        if type(s) == symGroup: s = self.addSymbol(sname)
        return s

    def moveGroup(self, name, group):
        """ move an existing symGroup in the SymbolTable"""
        if self.data.has_key(group.name): group = self.data.pop(group.name)
        owner,nam = self.__findName(name)
        fullname = "%s.%s" % (owner.name, nam)
        group.status = 'normal'
        if group.toplevel:
            fullname = nam
            owner = self.data
        owner[nam] = group
        group.setname(fullname)
        return group

    def addGroup(self,name,toplevel=False,status='normal'):
        """ add a Symbol Group, either to the current LocalGroup (the default)
        or at the toplevel.
        TODO: add toplevel / complex name??
                addGroup('_sys.newgroup')  works already
        """
        g = self.moveGroup(name,self.addTempGroup(toplevel=toplevel))
        g.status = status
        g.setname()
        return g

    def setVariable(self,name,value=None,**kws):
        "add a variable"
        return self.addSymbol(name,value=value,type=symTypes.variable,**kws)

    def setObject(self,name,value=None,**kws):
        "add a reference to a python object"
        return self.addSymbol(name,value=value,type=symTypes.object**kws)

    def setProcedure(self, name, code,desc=None, **kws):
        " add defined procedure "
        if desc is None: desc = name
        return self.addSymbol(name,value=name,type=symTypes.defpro,
                              mod=self.ModuleGroup,
                              code=code,desc=desc,**kws)
        
    def setDefVariable(self, name, code, desc,**kws):
        "add defined variable" 
        return self.addSymbol(name,value=None, type=symTypes.defvar,
                              code=code,desc=desc,**kws)

    def setFunction(self,name,func,ftype=symTypes.pyfunc,code=None,
                    desc=None,cmd_out=None):
        "add a function"
        if func is None: return None
        if ftype not in symTypes.Funcs:
            raise SymbolError, 'cannot add function %s with type  %s' % (name,ftype)

        fcn_kws = None
        try:
            if desc is None: desc = func.__doc__
            try:
                if (func.__name__.startswith('tdl') or 
                    'tdl' in inspect.getargspec(func)[0]):
                    fcn_kws = {'kws':{'tdl':self.tdl}}
            except TypeError:  # numpy ufuncs will raise a TypeError here...
                pass
        except:
            raise SymbolError, 'cannot add function %s ' % name

        if fcn_kws:
            return self.addSymbol(name,value=func,type=ftype,code=code,
                                  desc=desc,cmd_out=cmd_out,**fcn_kws)
        else:
            return self.addSymbol(name,value=func,type=ftype,code=code,
                                  desc=desc,cmd_out=cmd_out)


    def getSymbol(self,name):
        """ look up symbol name
        """
        names,rest = self.resolveName(name)
        sym = self.data[names[0]]
        for i in names[1:]:    sym = sym[i]
        # print 'getSymbol ', name, names,rest, self.data[names[0]], sym
        
        # 'rest' may contain attributes of an object
        for i in rest:
            sn = None
            try:
                sn = getattr(sym.value,i)
            except:
                try:
                    sn = getattr(sym,i)
                except:
                    pass
            if sn is None:
                if hasattr(sym,name):
                    msg = sym.name
                else:
                    msg = sym.__repr__()
                raise SymbolError, 'cannot find "%s" in "%s"' % (i,msg)
            sym = sn
        return sym
     
    def getSymbolValue(self,name):
        n = self.getSymbol(name)
        if n is not None: return n.value
        return None

    def getVariable(self,name):
        return self.getSymbol(name)

    def getFunction(self,name):
        "returns a function, assumin it exists"
        ret = None
        try:
            f = self.getSymbol(name)
            if isinstance(f,Symbol):
                if f.type in symTypes.Funcs: ret = f
            elif callable(f):
                ret = f
        except:
            pass
        return ret

    def getGroup(self,name):
        names,rest = self.resolveName(name)
        if rest != []:
            raise SymbolError, 'cannot get group %s' % name
        
        parent,child = self.data,names[0]
        for i in names[1:]:
            parent,child = parent[child],i
            if type(parent[child]) != symGroup:
                raise SymbolError, 'cannot get group %s' % name

        return (parent,child)
    
    def delGroup(self,name):
        names,rest = self.resolveName(name)
        if rest != []:
            raise SymbolError, 'cannot delete group %s' % name
        
        parent,child = self.data,names[0]
        for i in names[1:]:
            parent,child = parent[child],i
            if type(parent[child]) != symGroup:
                raise SymbolError, 'cannot delete group %s' % name

        if parent[child].status == 'nodelete':
            raise SymbolError, 'cannot delete group %s' % name            

        parent.pop(child)
        return
    
    def delVariable(self,name):
        names,rest = self.resolveName(name)
        if rest != []:
            raise SymbolError, 'cannot delete variable: "%s" -- python attribute??' % name
        
        parent,child = self.data,names[0]
        for i in names[1:]:
            parent,child = parent[child],i
        if isinstance(parent[child],Symbol):
            parent.pop(child)
        else:
            raise SymbolError, 'cannot delete variable %s -- group?? ' % name            
        return

    def delSymbol(self,name):
        names,rest = self.resolveName(name)
        if rest != []:
            raise SymbolError, 'cannot delete symbol: "%s" -- python attribute??' % name
        
        parent,child = self.data,names[0]
        for i in names[1:]:
            parent,child = parent[child],i

        parent.pop(child)
        return


    def hasSymbol(self,name):
        " returns whether a symbol exists or not"
        try:
            return (None != self.getSymbol(name))
        except SymbolError:
            return False


    def hasFunction(self,name):
        "returns whether a function exists"
        try:
            return (None != self.getFunction(name))
        except SymbolError:
            return False

    def hasGroup(self,name):
        "returns whether a Group exists"
        try:
            n = self.getSymbol(name)
            if type(n) == symGroup: return True
        except SymbolError:
            pass
        return False

    def initialize(self,libs=None):
        ll = init_loadLibs
        if libs is not None and type(libs) in (types.ListType,types.TupleType):
            ll = ll + list(libs)
        # 
        spath = self.getSymbolValue('_sys.path')
        if type(spath) == types.ListType:
            for p in spath: set_path(p)
        elif type(spath) == types.StringType:
            set_path(spath)

        # import all libraries
        for i in ll: self.import_lib(i)


            
    def import_lib(self,lib):
        " import or reload module given module name or object"
        if lib is None: return None

        mod, msg = None, None
        if type(lib) == types.StringType:
            try: 
                mod = __import__(lib)
                components = lib.split('.')
                for comp in components[1:]:
                    mod = getattr(mod, comp)
            except ImportError:
                msg = '    Error loading module %s:' % lib
        elif type(lib) == types.ModuleType:
            try:
                mod = reload(lib)
            except ImportError:
                msg = '    Error loading module %s:' % lib
        if mod is None:
            self.writer.write("    cannot load module %s !" % lib)
            if msg is not None: PrintShortExcept(msg)
            return None

        # mod is now a real module, not a string of the module name
        title = getattr(mod,'title',mod.__name__)
        self.writer.write("    loading %s ..." % title)
        self.writer.flush()
        if mod not in self.loaded_libs: self.loaded_libs.append(mod)
        
        try:
            for nam,val in getattr(mod,'_var_',{}).items():
                self.addSymbol(nam,val)
            for nam,val in getattr(mod,'_consts_',{}).items():
                self.addSymbol(nam,val,constant=True)
            for nam,val in getattr(mod,'_func_',{}).items():
                cmdOut = None
                func   = val
                if type(val) == types.TupleType:
                    func = val[0]
                    if len(val) > 1: cmdOut = val[1]
                x =self.setFunction(nam,func,cmd_out=cmdOut) 
            for nam in getattr(mod,'_scripts_',[]):
                try:
                    file_path = os.path.abspath(os.path.dirname(mod.__file__))
                    file_name = os.path.join(file_path,nam)
                    if os.path.exists(file_name) and os.path.isfile(file_name):
                        self.tdl.load_file(file_name)
                    else:
                        print "Warning: Cannot find lib script file: %s" % file_name
                except:
                    PrintExceptErr("Error loading script file '%s'"  % file_name)
            if self.tdl:
                for nam,val in getattr(mod,'_help_',{}).items():
                    self.tdl.help.add_topic(nam,val)
            import_msg = 'ok.'
        except ImportError:
            import_msg = 'import failed!'
        self.writer.write(" %s\n" % import_msg)

    def clearTempGroups(self):
        " clear all toplevel groups with delete status"
        k = self.data.keys()
        for i in k:
            if type(self.data[i]) == symGroup and self.data[i].status=='delete':
                self.data.pop(i)
                
    def addTempGroup(self,prefix='',nlen=8,**kw):
        """ add a random-ish group name, as for a procedures namespace """
        
        if prefix is None: prefix = ''
        grp  = symGroup(name='',status='delete',**kw)
        keys = self.data.keys()
        gname= "%s%s" % (prefix,self.randomName(n = nlen))
        ntry = 0
        while gname in keys:
            gname = "%s%s" % (prefix,self.randomName(n = nlen))                
            ntry = ntry + 1
            if (ntry>10000):
                raise SymbolError, 'cannot create a group name! %s' % prefix
            if (ntry%500 == 0):
                random.seed() ; time.sleep(0.001)
        self.data[gname] = grp
        grp.setname(gname)
        return grp
    
    def showTable(self,skip=None):
        print '======================'
        if skip is None:  skip = []
        def showGroup(x,prefix=None,indent=-1):
            px = '%s %s' % (' '*2*indent,prefix)
            for i in x.keys():
                if i not in skip:
                    s = i
                    if prefix is not None: s = "%s.%s" % (px,i)
                    print  s, ' ',  x[i].value
                    if type(x[i]) is symGroup:
                        showGroup(x[i], prefix=s,indent=indent+1)

        showGroup(self.data,prefix=None,indent=0)

    def listGroups(self):
        """collect full list of symbol groups """
        ret = []
        def showGroup(x,prefix=None):
            px = prefix
            if prefix is not None:
                ret.append(px)
            for i in x.keys():
                s = i
                if prefix is not None: s = "%s.%s" % (px,i)
                if type(x[i]) is symGroup:
                    showGroup(x[i], prefix=s)
        showGroup(self.data,prefix=None)
        return ret


    def saveTable(self,file):
        """ gather all pickle-able data types, save to file """
        store = copy.deepcopy(self.data)
        def removeFunctions(x):
            print x
        removeFunctions(store)
        

    def restoreTable(self,file):
        """ if group = None get all data, otherwise get all data in group """
        for sym in blob:
            g = sym.group
            n = sym.name
            if not self.hasGroup(g): self.addGroup(g)
            self.addVariable("%s.%s" % (g,n),value=deepcopy(sym.value),
                             type=sym.type,constant=sym.constant)
            
