
###########################################################
# Tdl Utilities Class
###########################################################


class wxTdlUtil:
    # make sure the below variables are set in the on_initialize method
    # then call init_tdl
    # self.shell = None
    # self.tdl = None

    # run this in the on_initialize method of the class
    def init_tdl(self):
        # get the tdl reference from the parent window 
        self.shell = self.getParent().get_shell()
        self.tdl   = self.shell.tdl

    def setValue(self,var_name,value):
        self.tdl.setVariable(var_name,value)

    def getValue(self,var_name):
        return self.tdl.getVariableValue(var_name)

    def getVariable(self,var_name):
        return self.tdl.getVariable(var_name)

    def listGroups(self):
        data = self.tdl.symbolTable.listGroups()
        print data.sort()
        return data

    def listDataGroup(self,grp):
        #return self.tdl.symbolTable.listDataGroup(grp)
        data = self.tdl.symbolTable.listGroups()
        print data.sort()
        return data

    def listAllData(self):
        #data_dict = self.tdl.symbolTable.listData()
        #lst = []
        #for grp in data_dict.keys():
        #    for node in data_dict[grp]:
        #        name = "%s.%s" % (grp,node)
        #        lst.append(name)
        data = self.tdl.symbolTable.listVariables()
        print data.sort()
        return data

    def execLine(self,line):
        #return self.tdl.execute(line)
        return self.tdl.eval(str(line))

    def str_to_list(self,s,conv=float):
        if s == None: return []
        if s == 'None': return []
        s = s.strip()
        if len(s) == 0: return []
        if s[0] == '=': s = s[1:]
        if s[0] == '[':
            s = s[1:]
            if s[-1] == ']':
                s = s[:-1]
        if len(s) == 0: return []
        if s.find(',') < 0:
            lst = s.split()
        else:
            lst = s.split(',')
        # we probably need to loop
        # and do explicit instead of use map
        #return map(conv, lst)
        v = []
        for l in lst:
            if conv == float:
                n = float(eval(l))
                v.append(n)
            elif conv == int:
                n = int(eval(l))
                v.append(n)
        return v

    def str_to_list_var(self,s,conv=float):
        if s == None: return []
        s = s.strip()
        if len(s) == 0: return []
        if s[0].isalpha() or s[0] == '_':
            v = self.getValue(s)
            if v == None: return []
            v2 = []
            for l in v:
                if conv == float:
                    n = float(l)
                    v2.append(n)
                elif conv == int:
                    n = int(l)
                    v2.append(n)
            return v2
        else:
            return self.str_to_list(s,conv=conv)

    def post_message(self,mess):
        """
        write a message to a component called Messages, or to screen
        """
        txt = "\n---\n%s\n" % mess
        try:
            #self.components.Messages.text = txt
            self.components.Messages.appendText(txt)
        except:
            print txt
        return
