import sys
from util import closure, randomName

MAX_GROUP_LEVEL = 16

def isGroup(g): return isinstance(g,Group)

class Group(object):
    """named container of other variables and subgroups: lightweight class"""
    def __init__(self):
        pass

    def __len__(self):
        return len(self.__dict__)

    def __repr__(self):
        return "<%s: %i items, id=%s>" % (self.__class__.__name__,
                                              len(self),hex(id(self)))
    def __setattr__(self,attr,val):
        """set group attributes."""
        self.__dict__[attr] = val

    def _members(self):
        "return sorted list of names of member"
        return sorted([i for i in self.__dict__.keys() if not i.startswith('_Group')])
    def _member_ids(self):
        "return dict of {id:attr} for member attributes"
        return dict([(id(i),i) for i in self.__dict__.keys()])


    def _subgroups(self):
        "return sorted list of names of members that are sub groups"
        return sorted([k for k,v in self.__dict__.items() if isGroup(v)])

class symbolTable(Group):
    top_group    = '_main'
    core_groups = ('_main','_sys','_builtin','_math')

    def __init__(self,tdl=None,writer=None,libs=None):
        Group.__init__(self)
        self.__writer = writer  or sys.stdout.write
        self.__inter  = {'local':None,
                         'module':None,
                         'path':None,
                         'localGroup':None,
                         'moduleGroup':None,
                         'searchGroups':None,
                         'createdGroup':False}
        self.__init_libs(libs=libs)
        
    def _set_internal(self,key,val):
        self.__inter[key] = val
        
    def __repr__(self):
        return "<%s: %i items, id=%s>" % (self.__class__.__name__,
                                          len(self),hex(id(self)))
    
    def __len__(self):
        return len(self.__dict__)-4

    def __init_libs(self,libs=None):
        sgroups = [self.top_group]

        for g in self.core_groups:
            if g == self.top_group or g in sgroups: continue
            setattr(self,g,  Group())
            sgroups.append(g)

        self._sys.path         = ['.']
        self._sys.searchGroups = sgroups
        self._sys.LocalGroup   = self.top_group
        self._sys.ModuleGroup  = self.top_group
        self._fix_searchGroups()
                
    def _import_libs(self,libnames=None):
        msg =  "IMPORT LIB NOT YET IMPlEMENTED: %s" %(libnames)
        self.__writer("%s\n" % msg)
        if libnames is not None:
            for i in libnames: self._import_lib(i)

    def _load_functions(self,funclist=None,group=None,parent=None,**kw):
        if group is None or funclist is None: return
        if isinstance(funclist,(list,tuple)) and \
           parent is not None:
            for name in funclist:
                setattr(group,name,getattr(parent,name))
        elif isinstance(funclist,dict):
            for name,fcn in funclist.items():
                setattr(group,name, closure(func=fcn,**kw))
    
    def _fix_searchGroups(self):
        """resolve list of groups to search for symbol names:

        The variable self._sys.searchGroups holds the list of group
        names for searching for symbol names.  A user can set this
        dynamically.  The names do need to be absolute (relative to
        _main).

        The variable self.__inter['searchGroups'] holds the list of 
        actual group objects resolved from this name.

        _sys.LocalGroup,_sys.ModuleGroup come first in the search list.
        """
        #
        # NOTE:  should look by id(var) not Group name!!
        #
        ##
        # check (and cache) whether searchGroups needs to be changed.
        if not (self._sys.LocalGroup   == self.__inter['local'] and
                self._sys.ModuleGroup  == self.__inter['module'] and
                self._sys.searchGroups == self.__inter['path']):
        
            print 'rebuilding search path'

            self.__inter['local']  = self._sys.LocalGroup
            self.__inter['module'] = self._sys.ModuleGroup 
            self.__inter['path']   = self._sys.searchGroups[:]


            # use current value of searchGroups to resolve group names
            if self.__inter['searchGroups'] is None:
                self.__inter['searchGroups'] = []
                
            old = self.__inter['searchGroups'][:]
            self.__inter['searchGroups'] = []        
            self._sys.searchGroups  = []

            for gname in (self.core_groups):
                if gname == self.top_group:
                    grp = self
                else:
                    grp = getattr(self,gname)
                if grp not in old: old.append(grp)

            loc_group = self._sys.LocalGroup or self.top_group
            mod_group = self._sys.ModuleGroup or self.top_group

            grps = [loc_group]
            if mod_group not in grps: grps.append(mod_group)
            for g in self._sys.searchGroups + list(self.core_groups):
                if g not in grps: grps.append(g)

            for name in grps:
                d = None
                if name == self.top_group:
                    d = self
                elif hasattr(self,name):
                    d = getattr(self,name)
                else:
                    for i in old:
                        if hasattr(i,name): d = getattr(i,name)
                if d is not None and d not in self.__inter['searchGroups']:
                    self.__inter['searchGroups'].append(d)
                    self._sys.searchGroups.append(name)

            for name,group in zip(self._sys.searchGroups,self.__inter['searchGroups']):
                if name == loc_group:  self.__inter['localGroup'] = group
                if name == mod_group:  self.__inter['moduleGroup'] = group

        return self.__inter

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
            names = g._members()
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
        self._fix_searchGroups()
        if name == self.top_group: return (self.top_group,self)

        parts = name.split('.')
        parts.reverse()
        top   = parts.pop()
        out   = None

        if top == self.top_group :
            out = self
            nam = [name]
        else:
            groupmap = zip(self.__inter['searchGroups'],
                           self._sys.searchGroups)

            for d,gnam in groupmap:
                if hasattr(d,top):
                    out = getattr(d,top)
                    parentname = gnam
                    nam = [top,parentname]
                    n  = 0
                    while parentname != self.top_group:
                        for d,gnam in groupmap:
                            if hasattr(d,parentname):
                                nam.append(gnam)
                                parentname = gnam
                        n = n+1
                        if n>MAX_GROUP_LEVEL:
                            raise LookupError, "cannot locate symbol '%s'" % name
                    break

        if out is None:
            raise LookupError, "cannot locate symbol '%s'" % name
        nam.reverse()
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
            nam.append(p)

        fullname =  '.'.join(nam)
        return fullname,out
        
    def getSymbol(self,s,create=False):
        "lookup and return a symbol by name"
        name,sym=self._lookup(s,create=create)
        return sym
      
    def getGroup(self,gname):
        "find group by name"
        name,sym=self._lookup(gname,create=False)
        if isGroup(sym):
            return sym
        else:
            raise LookupError, "symbol '%s' found, but not a group" % (name)

    def show_group(self,gname):
        print 'Show group ', gname, type(gname)
        if isGroup(gname): 
            grp = gname
            fullname = repr(grp)
        else:
            fullname,grp = self._lookup(gname,create=False)

        if not isGroup(grp):
            raise LookupError, "symbol '%s' found, but not a group" % (gname)

        title = fullname
        if title.startswith(self.top_group): title = title[6:]

        if grp == self: title = 'SymbolTable _main'

        mem = grp._members()
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
        

    def createGroup(self):
        return Group()
        
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
                setattr(grp,n,Group(n,status='normal'))

        setattr(grp,child,value)
        return getattr(grp,child)        

    def _parentOf(self,name):
        """return parent group, child name for an absolute symbol name (as from _lookup)
        that is, a pair suitable for hasattr, getattr, or delattr
        """
        n = name.split('.')
        if len(n) < 1 or name == self.top_group: return (self,None)
        child = n.pop()
        fname,sym=self._lookup('.'.join(n))
        return sym, child
    
    def delSymbol(self,name):
        name,sym=self._lookup(name,create=False)
        if isGroup(sym): 
            raise LookupError, "symbol '%s' is a group" % (name)
        parent,child = self._parentOf(name)
        if child is not None:  delattr(parent,child)

    def delGroup(self,name):
        name,sym=self._lookup(gname,create=False)
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

    print 'group1 members , subgroups: ', s.group1._members(), s.group1._subgroups()
# 
#     s.setSymbol('x',2.2)
#     s.setSymbol('ydef',DefinedVariable(expr='x*2.3',compiler=compiler))
#     print compiler.eval('ydef')
#     s.setSymbol('x',9)
#     print compiler.eval('ydef')
