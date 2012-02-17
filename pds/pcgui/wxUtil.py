"""
wX Utilities Class

Todo:
----
* make sure this independant of starting wxGUI....

"""
########################################################################

import sys

########################################################################
class wxUtil:
    # run this from the on_initialize method of the class
    def init_shell(self):
        self.shell = None
        # get the shell reference from the parent window
        # otherwise create a new one (ie for standalones)
        try:
            self.shell = self.GetParent().get_shell()
        except:
            from pds import shell
            self.shell = shell.Shell(stdin=self,stdout=self,GUI='WXAgg')

    ##############################################
    def set_data(self,var_name,value):
        self.shell.interp.set_data(var_name,value)

    ##############################################
    def get_data(self,var_name):
        return self.shell.interp.get_data(var_name)

    ##############################################
    def list_data(self,symbol=None,tunnel=False,_skip=True,instance=None):
        all = self.shell.interp.symbol_table.list_symbols(symbol=symbol,
                                                          tunnel=tunnel,
                                                          _skip=_skip,
                                                          instance=instance)
        dlist = all['var'] + all['ins']
        dlist.sort()
        return dlist

    ##############################################
    def exec_line(self,line):
        #return self.shell.interp.execute(str(line))
        return self.shell.exec_line(str(line))

    ##############################################
    def eval_line(self,line):
        return self.shell.interp.evaluate(str(line))

    ################################################
    def write(self, text):
        sys.stderr.write(text)
        
    ##############################################
    def post_message(self,mess):
        """
        write a message to a component called Messages, or to screen
        """
        txt = "\n--------------------\n"
        txt = txt + "%s\n--------------------\n" % mess
        try:
            #self.components.Messages.text = txt
            self.components.Messages.appendText(txt)
        except:
            #print txt
            self.write(txt)
        return
