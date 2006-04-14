# M. Newville Univ of Chicago (2005)
#
# -------------
# Modifications
# -------------
#
# 4-2-06 T2
# Modified symbolTable.initialize and symbolTable.import_lib
# to correctly handle strings as module names
#
# * 2-12-06 T2
# Added new methods for getting/putting to symbol table
# back to a single dictionary of symbols
#
##########################################################################

import types
import os
import random
import sys
import inspect
import re
from copy import deepcopy

from Num import Num
from Util import find_unquoted_char, split_delim, datalen
from Util import PrintExceptErr, PrintShortExcept

random.seed(0)

# should these be moved to a more common place?
DataTypes = ('variable','constant','defvar')
FuncTypes = ('pyfunc','defpro')
SymbolTypes = DataTypes + FuncTypes
SymbolTypeError = """SymbolTypeError: Valid type are:\n%s """ % str(SymbolTypes)
## NameError = """Invalid variable/function name"""

# Default data groups
#  search order will bw dataGroup,funcGroup,mainGroup,builtinGroup,mathGroup
builtinGroup   = '_builtin'
mainGroup    = '_main'
mathGroup    = '_math'
plotGroup    = '_plot'

isValidName = re.compile(r'[a-zA-Z_\$][a-zA-Z0-9_]*').match


class Symbol:
    """
    basic class for symbols for variables, functions, as stored in SymbolTable
    Symbols have a name, value, and type.  The type can be one of
    ** DATA
       constant   immutable variable (pi,e, etc)
       variable   regular variable
       defvar     defined variable (expression stack)
    ** FUNC
       pyfunc     python function
       defpro     defined procedure (sequence of statement code)
    """

    constant = False
    def __init__(self, name, value=None, desc=None, group=None,
                 help=None, type='variable',
                 code=None, constant = False,
                 as_cmd=True, cmd_out=None,
                 args=None, kws=None):
        self.name     = name
        self.constant = constant
        self.as_cmd   = as_cmd
        self.cmd_out  = cmd_out
        self.value    = value
        self.code     = code
        self.desc     = desc or ''
        self.help     = help or ''
        self.args     = args or []
        self.kws      = kws  or {}
        self.group    = group
        
        # check that type is valid
        if type not in SymbolTypes:
            raise SymbolTypeError

        self.type = type
        if self.type == 'constant': self.constant = True

    def getHelp(self):
        return self.help

    def getCode(self):
        return deepcopy(self.code)

    def __call__(self, *args,**kws):
        if self.type in ('pyfunc', 'pycmd'):
            x = {}
            x.update(self.kws) ; x.update(kws)
            val = self.value(*args,**x)
            return val
        else:
            return self.value

    def __cmdout__(self,val,**kws):
        if self.cmd_out:
            return self.cmd_out(val,**kws)
        elif val is None:
            return None
        return str(val)
        
    #def __setattr__(self, attr, val):
    #    if self.constant and  attr == 'value':
    #        print 'cannot set value for constant ', self.name
    #    else:
    #        self.__dict__[attr] = val

    def __repr__(self):
        name = "%s.%s" % (self.group,self.name)
        if self.type == 'variable':
            nelem = datalen(self.value)
            if nelem == 1:
                return "<Variable %s: type=scalar, value=%s>" % (name,repr(self.value))
            else: 
                if type(self.value) == Num.ArrayType:
                    return "<Variable %s: type=array, npts=%i, shape=%s>" % (name,len(self.value),self.value.shape)
                else:
                    t = str(type(self.value))[1:-1].replace('type','')
                    t.strip()
                    return "<Variable %s: type=%s, len=%i>" % (name,t,nelem)
                    
        elif self.type == 'constant':
            return "<Constant %s: value=%s>" % (name,repr(self.value))
        elif self.type == 'defvar':
            return "<Defined Variable %s ='%s'>" % (name,self.desc)
        elif self.type == 'pyfunc':
            return "<Function %s>" % (name)
        elif self.type == 'defpro':
            args = ','.join(self.args)
            for k,v in  self.kws.items(): args = "%s,%s=%s" % (args,k,str(v))
            return "<Procedure %s: args='%s'>" % (name,args)
        return "<Symbol %s: %s : %s>" % (name, self.type,repr(self.value))
    

class SymbolTable:
    """
    table of symbols and namespaces storing all functions and variables
    """
    def __init__(self,libs=None,writer=sys.stdout,tdl=None,**kws):

        self.tdl    = tdl
        self.writer = writer

        self.load_libs = []
        init_libs = []
        init_libs = ['TdlBuiltins','TdlNumLib','Plotter','IO']
        if libs is not None:
            init_libs.extend(libs)
        self.initialize(init_libs,clearAll=True)

    def initialize(self,libs=None,clearAll=False):
        if clearAll:
            self.sym    = {builtinGroup: {},
                           mathGroup: {},
                           plotGroup: {},                           
                           mainGroup: {}}
            self.dataGroup = mainGroup
            self.funcGroup = mainGroup
            self.addBuiltin('data_group',mainGroup)
            self.addBuiltin('func_group',mainGroup)
            self.setSearchGroups()
        if libs is not None:
            for lib in libs: self.import_lib(lib)

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
        if mod not in self.load_libs: self.load_libs.append(mod)
        
        try:
            for nam,val in getattr(mod,'_consts_',{}).items():
                self.addVariable(nam,val,constant=True)
            for nam,val in getattr(mod,'_var_',{}).items():
                self.addSymbol(nam,value=val,type='variable')
            for nam,val in getattr(mod,'_func_',{}).items():
                cmdOut = None
                asCmd  = True
                func   = val
                if type(val) == types.TupleType:
                    func = val[0]
                    if len(val) > 1: cmdOut = val[1]
                    if len(val) > 2: asCmd  = val[3]
                x =self.addFunction(nam,func,cmd_out=cmdOut,as_cmd=asCmd)
            if self.tdl:
                for nam,val in getattr(mod,'_help_',{}).items():
                    self.tdl.help.add_topic(nam,val)
            import_msg = 'ok.'
        except ImportError:
            import_msg = 'import failed!'                

        self.writer.write(" %s\n" % import_msg)

    ### Name/type/util functions
    def parseName(self,name,group=None,use_default=True):
        """
        split symbol name into (group,name) tuple.

        if name is fully qualified (that is, includes a prefix as 'group.name'),
        (group,name) are simply returned

        if the name is not fully qualified (no '.' in the name), the supplied group name
        is used, if available.

        if the group name is not provided and use_default=True, the current dataGroup name
        will be supplied as group.

        common cases:
           g,n = parseName('dat1.x')               => 'dat1', 'x'
           g,n = parseName('x')                    => self.dataGroup, 'x'
           g,n = parseName('x',group='dat2')       => 'dat2', 'x'
           g,n = parseName('x',use_default=False)  => None, 'x'
        """
        try:
            name.strip()
        except (AttributeError,NameError):
            raise NameError, 'cannot use symbol name %s ' % name        
        
        if '.' in name:
            try:
                group,name=name.split('.')
                group.strip()
                name.strip()
            except:
                raise NameError, 'cannot use symbol name %s ' % name

        # default group name for group == None
        # always defaults to current data group
        if group is None and use_default: group = self.dataGroup
        # print 'parseName ', group, name
        return (group, name)

    ### symbol manipulation functions
    def addSymbol(self,name,value=None,group=None,type='variable',code=None,desc=None,**kws):
        """
        add generic symbol to symbol table
        to specify which group the symbol goes to, you can either use
        name = group.name or  use name=name, group=group
        """
        # print 'Add Symbol ', name, value, code        
        group, name = self.parseName(name,group=group)
        if isValidName(group) and isValidName(name):
            if group not in self.sym.keys():self.sym[group]={}
            if name in self.sym[group].keys():
                if self.sym[group][name].constant:
                    #print "%s.%s is Constant type, cannot overwrite" % (group,name)
                    return (None,None)
            self.sym[group][name] = Symbol(name,value=value,type=type,code=code,
                                           desc=desc,group=group,**kws)
            return (group,name)
        else:
            return (None,None)

    def deleteSymbol(self,name,group=None,override=False):
        """
        delete generic symbol (py function or variable) in symbol table.
        to specify which group you can either use
        name = group.name or  use name=name, group=group
        """
        try:
            group, name = self.hasSymbolName(name,group=group)
        except NameError:
            return
        
        
        if group in (None, builtinGroup): return
        if self.sym[group][name].constant and not override:   return

        self.sym[group].pop(name)
        if len(self.sym[group]) == 0 and group not in (mainGroup,builtinGroup):
            self.deleteGroup(group)
            if self.dataGroup == group: self.dataGroup = mainGroup

    def getSymbol(self,name,groups=None,create=False):
        """
        return symbol (not just value!), creating if necessary
        name can be of form group.name or use just name and look in the list of supplied groups
        If type is in DataTypes default groups are [self.groupName,globalGroup,builtinGroup]
        If type is in FuncTypes default groups are [funcGroup]        
        creation puts symbol in the first group listed
        """
        # get group,name
        group, name = self.parseName(name,group=None,use_default=False)
        create_group = group

        # see if it exists, ie simple case with full name qualification:
        if group and self.sym.has_key(group):
            if self.sym[group].has_key(name):
                return self.sym[group][name]

        # if name was not fully qualified, then group=None, and
        # we need to search through groups
        if group is None:
            # search groups for name
            if groups is None:  groups = self.searchGroups
            create_group = groups[0]
            for group in groups:
                if self.sym[group].has_key(name):
                    return self.sym[group][name]

        if create:
            group, name = self.addSymbol(name,value=None,group=create_group,type='variable')
            return self.sym[group][name]
        return None

    def getSymbolValue(self,name,groups=None,default=None):
        sym = self.getSymbol(name,groups=groups,create=False)
        if sym:
            return sym.value
        else:
            return default


    #### Group manipulation and symbol checking
    def hasSymbol(self,name,group=None):
        "see if a symbol exists"
        try:
            group, name = self.parseName(name,group=group,use_default=False)
        except NameError:
            return False
        if self.hasGroup(group):
            return name in self.sym[group].keys()

        for group in self.searchGroups:
            if name in self.sym[group].keys(): return True
        return False

    def hasSymbolName(self, name, group=None):
        "  ret (group,name) or (None,None) if not exist "
        try:
            group, name = self.parseName(name,group=group,use_default=False)
        except NameError:
            return (None,None)
        if self.hasGroup(group):
            if name in self.sym[group].keys():   return (group,name)
        else:
            for grp in self.searchGroups:
                if name in self.sym[grp].keys(): return (grp,name)
        return (None,None)

    def SymbolType(name,group=None):
        try:
            group, name = self.hasSymbolName(name,group=group)
        except NameError:
            return (None,None,None)
        if group:
            ty = self.sym[group][name].type
            return (group, name, ty)
        else:
            return (None, None, None)

    def hasGroup(self, group):
        " does this group exist? "
        return (self.sym.has_key(group) and  group is not None)

    def addGroup(self, group):
        " add a Group "
        group = group.strip()
        if isValidName(group):
            if not self.sym.has_key(group): self.sym[group] = {}
        return group

    def deleteGroup(self, group):
        " delete a group and all its symbols (except default groups)"
        group = group.strip()
        #if not group in (builtinGroup,globalGroup) and self.sym.has_key(group):
        if not group in (builtinGroup,) and self.sym.has_key(group):
            self.sym.pop(group)
            if group == mainGroup:
                self.sym[mainGroup] = {}
        return None

    def getDataGroup(self,group=None,create=True):
        "return group name, or default (self.dataGroup), and makes sure the group exist"
        if group is None: group = self.dataGroup
        if create:        group = self.addGroup(group)
        return group

    def getFuncGroup(self,group=None,create=True):
        "return group name, or default (self.dataGroup), and makes sure the group exist"
        if group is None: group = self.funcGroup
        if create:        group = self.addGroup(group)
        return group

    def setSearchGroups(self):
        self.searchGroups = [self.dataGroup]
        for i in (self.funcGroup,mainGroup,builtinGroup,mathGroup,plotGroup):
            if i not in self.searchGroups: self.searchGroups.append(i)
        return self.searchGroups
    
    def getSearchGroups(self):
        return self.searchGroups
            
    def setDataGroup(self, group):
        " set the current active group "
        group = group.strip()
        if not self.hasGroup(group): self.sym[group] = {}
        self.dataGroup = group
        self.setSearchGroups()        
        self.setBuiltin('data_group',group)
        return group

    def setFuncGroup(self, group):
        " set the current active group "
        group = group.strip()
        if not self.hasGroup(group):  self.sym[group] = {}
        self.funcGroup = group
        self.setSearchGroups()        
        self.setBuiltin('func_group',group)
        return group

    def addRandomGroup(self,prefix='',nlen= 8):
        """
        add a quasi-randomly generated group name that is
        'guaranteed' to be unique.
        If prefix is provided, the group name begins with
        that string
        returns the group name or None.
        """
        rand = random.randrange
        def randomName(n=8):
            chars = 'abcdefghijklmnopqrstuvwxyz0123456789'
            return "_%s" % ''.join([chars[rand(0,len(chars))] for i in range(n)])
        if prefix is None: prefix = ''
        group = "%s%s" % (prefix,randomName(n = nlen))
        if self.hasGroup(group):
            ntry = 0
            while ntry<100:
                group = "%s%s" % (prefix,randomName(n = nlen))                
                if not self.hasGroup(group):break
                ntry = ntry + 1
                if (ntry % 20 == 0): nlen = nlen+1
        if self.hasGroup(group):
            return None
        self.sym[group] = {}
        return group
    
    # generic add/delete/get data
    def hasData(self,name,group=None):
        "add data"
        try:
            group, name = self.hasSymbolName(name,group=group)
        except NameError:
            return False
        if group:
            if self.sym[group][data].type in DataTypes:
                return True
        else:
            return False

    def addData(self,name,value=None,code=None,type='variable'):
        "add data"
        if type not in DataTypes: raise SymbolTypeError
        return self.addSymbol(name,value=value,group=self.groupName,code=code,type=type)
        
    def deleteData(self,name):
        "delete a data entry"
        (group, name, type) = self.SymbolType(name)
        if group and name:
            if type in DataTypes:
                return self.deleteSymbol(name)

    def getData(self,name):
        "get data symbol"
        sym = self.getSymbol(name,groups=None,create=False)
        if sym is None: return None
        if sym.type in DataTypes:
            return sym
        else:
            return None

    def getAllData(self,group=None):
        """ if group = None get all data, otherwise get all data in group """
        data = []
        if group is None:
            for group in self.sym.keys():
                for name in self.sym[group].keys():
                    if self.sym[group][name].type in DataTypes:
                        data.append(self.sym[group][name])
        else:
            if self.hasGroup(group):
                for name in self.sym[group].keys():
                    if self.sym[group][name].type in DataTypes:
                        data.append(self.sym[group][name])
        return data

    def clearAllData(self,group=None):
        " delete data "
        # work on how this should operate...
        doomed = []
        if group is None:
            doomed = self.sym.keys()
        elif self.hasGroup(group):
            doomed = [group]
        for group in doomed:
            for name in self.sym[group].keys:
                if self.sym[group][name].type in DataTypes:
                    self.sym[group].pop(name)
            if len(self.sym[group]) == 0:
                self.deleteGroup(group)

    # generic add/delete/get func
    def hasFunc(self,name,group=None):
        try:
            group, name = self.hasSymbolName(name,group=group)
        except NameError:
            return False
        if group is None: return False
        return self.sym[group][name].type in FuncTypes

    def hasFuncAsCmd(self,name,group=None):
        ' ever used any more??'
        try:
            group, name = self.hasSymbolName(name,group=group)
        except NameError:
            return False
        if group is None: return False
        if self.sym[group][name].type in FuncTypes:
            return self.sym[group][name].as_cmd        

    def addFunction(self,name,func,ftype='pyfunc',code=None,
                    desc=None,as_cmd=True,cmd_out=None):
        "add function"
        if not func: return (None,None)
        if ftype not in FuncTypes:
            raise SymbolTypeError, 'cannot add function %s ' % name

        fcn_kws = None
        try:
            if desc is None: desc = func.__doc__
            # look for 'tdl' argument --
            #    if present, add kw-arg to pass to function
            try:
                if 'tdl' in inspect.getargspec(func)[0]:
                    fcn_kws = {'kws':{'tdl':self.tdl}}
            except TypeError:  # numpy ufuncs will raise a TypeError here...
                pass
        except:
            raise SymbolTypeError, 'cannot add function %s ' % name

        if fcn_kws:
            return self.addSymbol(name,value=func,group=self.funcGroup,
                                  type=ftype,code=code,desc=desc,
                                  as_cmd=as_cmd,cmd_out=cmd_out,**fcn_kws)
        else:
            return self.addSymbol(name,value=func,group=self.funcGroup,
                                  type=ftype,code=code,desc=desc,
                                  as_cmd=as_cmd,cmd_out=cmd_out)
            
    def deleteFunc(self,name):
        "delete a func entry"
        (group,name,type) = self.SymbolType(name)
        if group and name:
            if type in FuncTypes:
                return self.deleteSymbol(name)

    def getFunc(self,name):
        "get func symbol"
        sym = self.getSymbol(name,groups=None,create=False)
        if sym is None: return None
        if sym.type in FuncTypes:
            return sym
        else:
            return None

    def getAllFunc(self,group=None):
        """ if group = None get all data, otherwise get all data in group """
        func = {}
        if group is None:
            for group in self.sym.keys():
                for name in self.sym[group].keys:
                    if self.sym[group][name].type in FuncTypes:
                        x = "%s.%s" % group,name
                        func.update({x:sym[group][name]}) 
        else:
            if self.hasGroup(group):
                for name in self.sym[group].keys:
                    if self.sym[group][name].type in FuncTypes:
                        x = "%s.%s" % group,name
                        func.update({x:sym[group][name]})
        return func

    def clearAllFunc(self,group=None):
        " delete data "
        doomed = []
        if group is None:
            doomed = self.sym.keys()
        elif self.hasGroup(group):
            doomed = [group]

        for group in doomed:
            for name in self.sym[group].keys():
                if self.sym[group][name].type in FuncTypes:
                    self.sym[group].pop(name)
            if len(self.sym[group]) == 0:
                self.deleteGroup(group)

    # variables (and const)
    def addVariable(self,name,value=None,constant=False):
        " add variable (or const)"
        if constant == True:
            type='constant'
        else:
            type='variable'
        group,name = self.addSymbol(name,value=value,type=type)
        return (group,name)

    def getVariable(self,name,create=False):
        """ get variable (or const) """
        sym = self.getSymbol(name,create=create)
        if sym is not None:
            if sym.type not in ('variable','constant','defvar'):
                sym = None
        return sym

    def getVariableCurrentGroup(self,name):
        """
        find variable (or const) in current group (or in group.name) or create it,
        WITHOUT looking to _builtin or _main or any other groups
        this should be used to lookup symbols for assignments
        """
        sym = self.getSymbol(name,groups=(self.dataGroup,),create=True)
        if sym is not None:
            if sym.type not in ('variable','constant','defvar'):
                sym = None
        return sym

    def setVariable(self,name,value,group=None):
        return self.addSymbol(name,value=value,type='variable',group=group)

    # builtins
    def addBuiltin(self,name,value,desc=None):
        " add a symbol to the _builtin group, and make it a 'constant'"
        if name.find('.') > -1:
            raise NameError, name
        return self.addSymbol(name,value=value,group=builtinGroup,type='constant')

    def setBuiltin(self,name,value,desc=None):
        " set a Builtin value "
        if name.find('.') > -1: raise NameError
        if not self.hasGroup(builtinGroup): self.sym[builtinGroup] = {}
        if not self.sym[builtinGroup].has_key(name):
            self.addBuiltin(name,value,desc=desc)
        else:
            self.sym[builtinGroup][name].value    = value
            if desc: 
                self.sym[builtinGroup][name].desc = desc
            self.sym[builtinGroup][name].constant = True
            self.sym[builtinGroup][name].type = 'constant'
        # if setting group name, make sure it exists in data
        if name == 'group':
            if not self.sym.has_key(value):  self.sym[value] = {}
            self.groupName = value
        return    

    # def variables
    def setDefVariable(self, name, code, desc,**kws):
        # print 'SetDEF VAR  ', name, code
        return self.addSymbol(name,value=None,type='defvar',code=code,desc=desc)
    
    # defined procedures    
    def addDefPro(self, name, code,desc=None, **kws):
        " add defined procedure "
        if desc is None: desc = name
        return self.addSymbol(name,value=name,group=self.funcGroup,
                              type='defpro',code=code,desc=desc,**kws)

    # functions to list stuff
    def listGroups(self):
        s = self.sym.keys()
        s.sort()
        return s

    def listSymbols(self):
        " return a dict of symbol names"
        all = {}
        groups = self.listGroups()
        for g in groups:
            sym = self.listGroupsSymbols(g)
            all.update({g:sym})
        return all

    def listGroupSymbols(self,group):
        if self.hasGroup(group):
            s = self.sym[group].keys()
            s.sort()
            return s
        else:
            return []

    def listData(self):
        "return a dict of data names as {grp:[name,name,...]}"
        all = {}
        glist = self.listGroups()
        for g in glist:
            s = self.listGroupSymbols(g)
            d = []
            for name in s:
                if self.sym[g][name].type in DataTypes:
                    d.append(name)
            all.update({g:d}) 
        return all

    def listFunc(self):
        "return a dict of func names as {grp:[name,name,...]}"
        all = {}
        glist = self.listGroups()
        for g in glist:
            s = self.listGroupSymbols(g)
            fcn = []
            for name in s:
                if self.sym[g][name].type in FuncTypes:
                    fcn.append(name)
            all.update({g:fcn})
        return all

    def listDataGroup(self,group=None):
        "return a list of Data names in a group"
        if group is None: group = self.getDataGroup(group)
        all = self.listData()
        if group in all.keys():
            return all[group]
        else:
            return []
    
    def listFuncGroup(self,group=None):
        "return a list of Func names in a group"
        if group is None: group = self.getFuncGroup(group)
        all = self.listFunc()
        if group in all.keys():
            return all[group]
        else:
            return []

