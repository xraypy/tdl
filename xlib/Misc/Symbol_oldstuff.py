
u = '''
class Old:
    
    def addGroup(self, gname):
        " add a Group "
        gname= gname.strip()

        if isValidName(group):
            if not self.groups.has_key(group): self.groups[group] = symGroup(name=group)
        return group

    ### symbol manipulation functions
    def addSymbol(self,name,value=None,group=None,type=symTypes.variable,
                  code=None,desc=None,constant=False,**kws):
        """add symbol to symbol table
        to specify which group the symbol goes to, you can either use
        name = group.name or  use name=name, group=group
        """

    def addRandomGroup(self,prefix='',nlen= 8):
        """
        add a quasi-randomly generated group name that is 'guaranteed' to be unique.
        If prefix is provided, the group name begins with that string
        returns the group name.
        """
        if prefix is None: prefix = ''
        group = "%s%s" % (prefix,randomName(n = nlen))
        if self.hasGroup(group):
            ntry  = 0
            while self.hasGroup(group):
                group = "%s%s" % (prefix,randomName(n = nlen))                
                ntry = ntry + 1
                if (ntry%10000 == 0):  
                    if (ntry>100000):
                        raise SymbolError, 'cannot create a group name! %s' % prefix
                    random.seed()
        self.groups[group] = symGroup(name=group)
        return group

    def delSymbol(self,name,group=None,override=False):
        """
        delete generic symbol (py function or variable) in symbol table.
        to specify which group you can either use
        name = group.name or  use name=name, group=group
        """
        group, name = self.hasSymbolName(name,group=group)
        if group in (None, builtinGroup): return
        if self.groups[group][name].constant and not override:   return

        self.groups[group].pop(name)
        if len(self.groups[group]) == 0 and group not in requiredGroups:
            self.delGroup(group)
            if self.dataGroup == group: self.dataGroup = mainGroup

    def getSymbol(self,name,groups=None,create=False):
        """
        return symbol (not just value!), creating if necessary
        name can be of form group.name or use just name and look in the list of supplied groups
        If type is in Datas default groups are [self.groupName,globalGroup,builtinGroup]
        If type is in Func default groups are [funcGroup]        
        creation puts symbol in the first group listed
        """
        # get group,name
        # print 'getSymbol ', name
        group, name = split_name(name,group=None)
        create_group = group

        # see if it exists, ie simple case with full name qualification:
        if group and self.groups.has_key(group):
            if self.groups[group].has_key(name):
                return self.groups[group][name]

        # if name was not fully qualified, then group=None, and
        # we need to search through groups
        if group is None:
            # search groups for name
            if groups is None:  groups = self.searchGroups
            create_group = groups[0]
            for group in groups:
                if self.groups[group].has_key(name):
                    return self.groups[group][name]

        if create:
            group, name = self.addSymbol(name,value=None,group=create_group,type=symTypes.variable)
            return self.groups[group][name]
        return None

    def getSymbolValue(self,name,groups=None,default=None):
        sym = self.getSymbol(name,groups=groups,create=False)
        if sym:
            return sym.value
        else:
            return default

    #### Group manipulation and symbol checking
    def hasSymbol(self,name,group=None):
        " returns whether a symbol exists or not"
        return (None,None) != self.hasSymbolName(name,group=group)

    def hasSymbolName(self, name, group=None):
        " sees if a symbol exists, returning (group,name) if it does exist, or (None,None) if not."
        try:
            group, name = split_name(name,group=group)
        except NameError:
            return (None,None)
        if self.hasGroup(group):
            if name in self.groups[group].keys():   return (group,name)
        else:
            for g in self.searchGroups:
                if name in self.groups[grp].keys(): return (grp,name)
        return (None,None)


    def hasGroup(self, group):
        " does this group exist? "
        return (self.groups.has_key(group) and  group is not None)

    def delGroup(self, group):
        " delete a group and all its symbols (except default groups)"
        group = group.strip()
        #if not group in (builtinGroup,globalGroup) and self.groups.has_key(group):
        #if self.groups.has_key(group) and group not in requiredGroups:
        #    self.groups.pop(group)
        # should we allow deleting required groups and just re-instate?
        # ie except for immutables??
        if group in immutableGroups:
            print "'%s' can not be deleted" % group
            return None
        elif self.groups.has_key(group):
            self.groups.pop(group)
            if group in requiredGroups:
                self.addGroup(group)
        else:
            print "Cant find group '%s'" % group
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
        if self.funcGroup != self.dataGroup: self.searchGroups.append(self.funcGroup)
        for i in requiredGroups:
            if i not in self.searchGroups: self.searchGroups.append(i)
        return self.searchGroups
    
    def getSearchGroups(self):
        return self.searchGroups
            
    def setDataGroup(self, group):
        " set the current active group "
        group = group.strip()
        if not self.hasGroup(group): self.groups[group] = symGroup(name=group)
        self.dataGroup = group
        self.setSearchGroups()        
        self.setVariable('data_group',value=group,group=builtinGroup)
        return group

    def setFuncGroup(self, group):
        " set the current active group "
        group = group.strip()
        if not self.hasGroup(group):  self.groups[group] = symGroup(name=group)
        self.funcGroup = group
        self.setSearchGroups()        
        self.setVariable('func_group',value=group,group=builtinGroup)
        return group


    
    def getAllData(self,group=None):
        """ if group = None get all data, otherwise get all data in group """
        data = []
        grouplist = [group]
        if group is None: grouplist = self.groups.keys()
        for group in grouplist:
            for name in self.groups[group].keys():
                if self.groups[group][name].type in symTypes.Data:
                    data.append(self.groups[group][name])
        return data

    def restoreData(self,blob):
        """ if group = None get all data, otherwise get all data in group """
        for sym in blob:
            g = sym.group
            n = sym.name
            if not self.hasGroup(g): self.addGroup(g)
            self.addVariable("%s.%s" % (g,n),value=deepcopy(sym.value),
                             type=sym.type,constant=sym.constant)
            
    def clearAllData(self,group=None):
        " delete data "
        # work on how this should operate...
        doomed = []
        if group is None:
            doomed = self.groups.keys()
        elif self.hasGroup(group):
            doomed = [group]
        for group in doomed:
            for name in self.groups[group].keys:
                if self.groups[group][name].type in symTypes.Data:
                    self.groups[group].pop(name)
            if len(self.groups[group]) == 0:
                self.delGroup(group)

    # generic add/delete/get func
    def hasFunc(self,name,group=None):
        try:
            group, name = self.hasSymbolName(name,group=group)
        except NameError:
            return False
        if group is None: return False
        return self.groups[group][name].type in symTypes.Funcs

            
    def delFunc(self,name):
        "delete a function symbol"
        sym = self.getSymbol(name)
        if sym is not None:
            if sym.type in symTypes.Funcs:
                return self.delSymbol(name)

    def getFunc(self,name):
        "get a function symbol"
        sym = self.getSymbol(name,groups=None,create=False)
        if sym is None: return None
        if sym.type in symTypes.Funcs:
            return sym
        else:
            return None

    def getVariable(self,name,create=False):
        """ get variable (or const) """
        sym = self.getSymbol(name,create=create)
        if sym is not None: 
            if sym.type not in symTypes.Data:   sym = None
        return sym


    def getVariableCurrentGroup(self,name):
        """
        find variable (or const) in current group (or in group.name) or create it,
        WITHOUT looking to _builtin or _main or any other groups
        this should be used to lookup symbols for assignments
        """
        sym = self.getSymbol(name,groups=(self.dataGroup,),create=True)
        if sym is not None:
            if sym.type not in symTypes.Data:
                sym = None
        return sym

    def setVariable(self,name,value,type=symTypes.variable,group=None,constant=False):
        return self.addSymbol(name,value=value,type=type,
                              group=group,constant=constant)

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
        s = self.groups.keys()
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
            s = self.groups[group].keys()
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
                if self.groups[g][name].type in symTypes.Data:
                    d.append(name)
            all.update({g:d}) 
        return all

    def listFunctions(self):
        "return a dict of func names as {grp:[name,name,...]}"
        all = {}
        glist = self.listGroups()
        for g in glist:
            s = self.listGroupSymbols(g)
            fcn = []
            for name in s:
                if self.groups[g][name].type in symTypes.Funcs:
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

    # obsolete functions:
    def hasData(self,name,group=None):
        PrintShortExcept('symbolTable.hasData is obsolete')        

    def addData(self,name,value=None,code=None,**kw):
        PrintShortExcept('symbolTable.addData is obsolete')                
        
    def delData(self,name):
        PrintShortExcept('symbolTable.delData is obsolete')

    def getData(self,name):
        PrintShortExcept('symbolTable.getData is obsolete')        

    def SymbolType(name,group=None):
        PrintShortExcept('symbolTable.SymbolType is obsolete')
'''
