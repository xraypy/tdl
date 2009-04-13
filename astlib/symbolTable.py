from __future__ import print_function
import sys
from util import closure, randomName

MAX_GROUP_LEVEL = 16

def isGroup(g): return isinstance(g,Group)

class Group(object):
    """container for variables, subgroups, and modules:
    a lightweight object, with 
    """
    def __init__(self,_name=None,_file=None,_status=None):
        self.__name = _name
        self.__file = _file
        self.__status = _status
        
    def __len__(self):  return len(dir(self))

    def __repr__(self):
        stats = "%i items, id=%s" % (len(self),hex(id(self)))
        front = "Group"
        if self.__name is not None: front = "Group '%s'" % self.__name
        return "<%s: %s>" % (front,stats)

    def __setattr__(self,attr,val):
        """set group attributes."""
        self.__dict__[attr] = val

    def __dir__(self):
        "return sorted list of names of member"
        return sorted([i for i in self.__dict__.keys() if not i.startswith('_Group__')])
   
    def _subgroups(self):
        "return sorted list of names of members that are sub groups"
        return sorted([k for k,v in self.__dict__.items() if isGroup(v)])

    def _name(self):
        return self.__name

    def _file(self):
        return self.__file

    def _status(self):
        return self.__status


class invalidName:
    # used to create a value that will NEVER be a useful symbol.
    # symbolTable._lookup() uses this to check for invalid names
    pass

class symbolTable(Group):
    top_group   = '_main'
    core_groups = ('_main','_sys','_builtin','_math')
    mod_group   = '_modules'
    __invalid_name = invalidName()

    def __init__(self,tdl=None,writer=None,libs=None):
        Group.__init__(self,_name=self.top_group)

        self.__writer = writer  or sys.stdout.write
        self.__cache  = {'localGroup':None, 'moduleGroup':None,'searchGroups':None}

        setattr(self,self.top_group, self)
        setattr(self,self.mod_group, Group(_name=self.mod_group))
        for g in self.core_groups:
            if not hasattr(self,g):
                setattr(self,g,  Group(_name=g))

        self._sys.path         = ['.']

        self._sys.searchGroups = list(self.core_groups)
        self._sys.search       = getattr(self,self.mod_group)

        self._sys.LocalGroup   = self
        self._sys.ModuleGroup  = self
        self._sys.pymodules    = sys.modules
        self._sys.modules      = {}
        
        for i in dir(self):
            print("%s = %s" % (i, repr(getattr(self,i))))
        self._fix_searchGroups()
    
    def __init_libs(self,libs=None):
        pass
        
    def import_module(self,name, asname = None, fromlist = None):
       msg =  "import %s"  % name
       if asname is not None: msg =  "%s as %s" % (ms,asname)
       
       self.__writer("%s\n"% msg)
       # first look for "name.tdl"
       isTDLmodule = False

       if not isTDLmodule:
           __import__(name)


    def _load_functions(self,funclist=None,group=None,parent=None,**kw):
        if group is None or funclist is None: return
        if isinstance(funclist,(list,tuple)) and \
           parent is not None:
            for name in funclist:
                setattr(group,name,getattr(parent,name))
        elif isinstance(funclist,dict):
            for name,fcn in funclist.items():
                setattr(group,name, closure(func=fcn,**kw))
    
    def _set_local_mod(self,groups):
        self._sys.LocalGroup, self._sys.ModuleGroup  = groups
        self._fix_searchGroups()
        
    def _fix_searchGroups(self):
        """resolve list of groups to search for symbol names:

        The variable self._sys.searchGroups holds the list of group
        names for searching for symbol names.  A user can set this
        dynamically.  The names need to be absolute (relative to
        _main).

        The variable self.__cache['searchGroups'] holds the list of 
        actual group objects resolved from this name.

        _sys.LocalGroup,_sys.ModuleGroup come first in the search list,
        followed by any search path associated with that module (from
        imports for that module)
        """
        ##
        # check (and cache) whether searchGroups needs to be changed.
        sys = self._sys
        cache = self.__cache

        if not (sys.LocalGroup   == cache['localGroup'] and
                sys.ModuleGroup  == cache['moduleGroup'] and
                sys.searchGroups == cache['searchGroups']):

            if sys.ModuleGroup is None: sys.ModuleGroup = self.top_group
            if sys.LocalGroup is None: sys.LocalGroup = self.ModuleGroup

            cache['localGroup']  = sys.LocalGroup 
            cache['moduleGroup'] = sys.ModuleGroup

            cache['searchGroups'] = [sys.LocalGroup]
            if sys.ModuleGroup != sys.LocalGroup: cache['searchGroups'].append(sys.ModuleGroup)

            for group in self._modules._subgroups():
                if group not in cache['searchGroups']:
                    cache['searchGroups'].append(group)

            for gname in self.core_groups:
                group = getattr(self,gname)
                if group not in cache['searchGroups']:
                    cache['searchGroups'].append(group)

            sys.searchGroups = cache['searchGroups']
        
        # print("end of fix: ", cache)
        return cache

    def _get_localgroup(self):   return self._fix_searchGroups()['localGroup']
    def _get_modulegroup(self):  return self._fix_searchGroups()['moduleGroup']    
    def _get_searchGroups(self): return self._fix_searchGroups()['searchGroups']

    def list_groups(self,group=None):
        if group in (self.top_group,None):
            g = self
            group = 'SymbolTable'
        elif hasattr(self,group):
            g = getattr(self,group)
        else:
            g = None
            msg = '%s not found' % group
            
        if isGroup(g):
            names = dir(g)
            o = ['== %s ==' % group]
            for i in names:
                if not (i.startswith('_Group__') or
                        i.startswith('_symbolTable__')):
                    o.append('  %s: %s' % (i,repr(getattr(g,i))))
            msg = '\n'.join(o)
        else:
            msg = '%s is not a Subgroup' % group
        self.__writer("%s\n" % msg)

    def has_symbol(self,s,group=None):
        "return whether there is a toplevel symbol with the give name"
        return hasattr(self,s)

    def has_group(self,s,group=None):
        "return whether there is a toplevel group with given name"
        return hasattr(self,s) and isGroup(getattr(self,s))
        
    def _lookup(self,name=None,create=False):
        """looks up symbol in search path
        returns (absolute_symbol_name,object) for symbol"""
        # print(" ... lookup name = %s " % name)

        self._fix_searchGroups()
        if name == self.top_group: return (self.top_group,self)

        parts = name.split('.')
        parts.reverse()
        top   = parts.pop()
        out   = self.__invalid_name
        if top == self.top_group:
            out = self
        else:
            grouplist = self.__cache['searchGroups']
            for grp in grouplist:
                if hasattr(grp,top):
                    out = getattr(grp,top)

        if out is self.__invalid_name:
            raise LookupError, "cannot locate symbol(1) '%s'" % name

        while parts:
            p = parts.pop()
            if hasattr(out,p):
                out = getattr(out,p)
            elif create: 
                val = None
                if len(parts) > 0: val = Group()
                setattr(out,p,val)
                out = getattr(out,p)
            else:
                raise LookupError, "cannot locate member '%s' of '%s'" % (p,out)

        return out
        
    def getSymbol(self,s,create=False):
        "lookup and return a symbol by name"
        return self._lookup(s,create=create)
      
    def getGroup(self,gname):
        "find group by name"
        sym=self._lookup(gname,create=False)
        if isGroup(sym):
            return sym
        else:
            raise LookupError, "symbol '%s' found, but not a group" % (gname)

    def show_group(self,gname):
        # print('Show group ', gname, type(gname))
        if isGroup(gname): 
            grp = gname
            title = grp._name() or 'Group'
        else:
            grp = self._lookup(gname,create=False)
            title = gname

        if not isGroup(grp):
            raise LookupError, "symbol '%s' found, but not a group" % (gname)

        if title.startswith(self.top_group): title = title[6:]

        if grp == self: title = 'SymbolTable _main'

        mem = dir(grp)
        o = ['== %s : %i symbols ==' % (title,len(mem))]
        for i in mem:
            if not (i.startswith('_Group__') or
                    i.startswith('_symbolTable__')):
                o.append('  %s: %s' % (i,repr(getattr(grp,i))))
        msg = '\n'.join(o)
        self.__writer("%s\n" % msg)


    def placeGroup(gname,group=None,parent=None):
        if parent is None: parent = selg._get_localgroup()
        if isinstance(parent,(str,unicode)):
            parent = self.getGroup(parent)
        if group is None:
            group = self.getGroup(gname)

        if isGroup(group):
            setattr(parent,gname,group)
        

    def createGroup(self,**kw):    return Group(**kw)
        
    def setSymbol(self,name,value=None,group=None):
        grp = self._get_localgroup()
        if group is not None:  grp = self.getGroup(group)
        names= name.split('.')
        child = names.pop()
        for n in names:
            if hasattr(grp,n):
                grp = getattr(grp,n)
                if not isGroup(grp):
                    raise ValueError, "cannot create subgroup of non-group '%s'" % grp
            else:
                setattr(grp,n,Group(_name=n)) ##  ,status='normal'))

        setattr(grp,child,value)
        return getattr(grp,child)        

    def _parentOf(self,name):
        """return parent group, child name for an absolute symbol name
        (as from _lookup) that is, a pair suitable for hasattr,
        getattr, or delattr 
        """
        n = name.split('.')
        if len(n) < 1 or name == self.top_group: return (self,None)
        child = n.pop()
        sym=self._lookup('.'.join(n))
        return sym, child
    
    def delSymbol(self,name):
        sym=self._lookup(name,create=False)
        if isGroup(sym): 
            raise LookupError, "symbol '%s' is a group" % (name)
        parent,child = self._parentOf(name)
        if child is not None:  delattr(parent,child)

    def delGroup(self,name):
        sym=self._lookup(gname,create=False)
        if not isGroup(sym): 
            raise LookupError, "symbol '%s' found, but not a group" % (name)
        if sym._Group__status == 'nodelete':
            self.__writer("cannot delete group '%s'\n" % name)
        else:
            parent,child = self._parentOf(name)
            if child is not None:  delattr(parent,child)
            
if __name__ == '__main__':
    s = symbolTable()
    from compiler import Compiler, DefinedVariable
    compiler = Compiler(symtable=s)

    s.group1 = Group()
    s.group2 = Group()

    s.show_group('_sys')
    s.group1.x = 12.0
    s.group1.g1 = Group()

    s.show_group('group1')
    s.group1.g1.title = 'a string here'
    s.group1.g1.x = 99120.102
    s.group1.g1.e = 8980.0
    
    s.show_group('group1.g1')
    s.list_groups()

    print('group1 members , subgroups: ', dir(s.group1), s.group1._subgroups())
# 
#     s.setSymbol('x',2.2)
#     s.setSymbol('ydef',DefinedVariable(expr='x*2.3',compiler=compiler))
#     print( compiler.eval('ydef'))
#     s.setSymbol('x',9)
#     print( compiler.eval('ydef'))
