
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
        return self.tdl.symbolTable.listGroups()

    def listDataGroup(self,grp):
        return self.tdl.symbolTable.listDataGroup(grp)

    def listAllData(self):
        data_dict = self.tdl.symbolTable.listData()
        lst = []
        for grp in data_dict.keys():
            for node in data_dict[grp]:
                name = "%s.%s" % (grp,node)
                lst.append(name)
        lst.sort()
        return lst

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
        return map(conv, lst)

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
