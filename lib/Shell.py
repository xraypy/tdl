#!/usr/bin/env python
############################################################################
#
# -------------
# Modifications
# -------------
# 4-8-06 T2
# small tweak to setting stdin/stdout so tdl_wxGUI will work.
# small change to do_show_symbols, checks if key is a group.
#
############################################################################

from Num import Num

from Eval   import Evaluator
from Symbol import SymbolTable
from Util   import ParseException, EvalException, show_list, split_arg_str, show_more
from Help   import Help
import version
import cmd
import os
import sys
import types
import getopt
import traceback


class shell(cmd.Cmd):
    banner = "\n  Tiny Data Language %s  M. Newville, T. Trainor (2006)"
    intro  = "\n  type 'help' to get started"
    ps1    = "tdl> "
    ps2    = "...> "
    max_save_lines = 500

    def __init__(self, completekey='tab', libs=None, debug=False,
                 stdin=None, stdout=None, intro=None, GUI='TkAgg'):

        cmd.Cmd.__init__(self,completekey='tab')

        if stdin is not None:
            sys.stdin = self.stdin = stdin
            self.use_rawinput = False
        else:
            self.stdin = sys.stdin
            self.use_rawinput = True
        if stdout is not None:
            sys.stdout = self.stdout = stdout
            #sys.stderr = stdout
        else:
            self.stdout = sys.stdout
        #self.stdin  = stdin  or sys.stdin
        #self.stdout = stdout or sys.stdout

        print self.banner % version.version

        self.tdl       = Evaluator(interactive= True, debug= debug,
                                  output = self.stdout, libs= libs, GUI=GUI)


        self.help      = self.tdl.help.get_help
        self.help_topics = self.tdl.help.list_topics().split()
        self.help_topics.sort()

        self.prompt    = self.ps1
        self.tdl.prompt= self.ps2
        self._status   = True
        self.debug     = debug
        
    def emptyline(self):
        pass
             
    def show_function(self,args):
        lout = ''
        if len(args)==0:
            f = self.tdl.symbolTable.listFunc()
            grps = f.keys()
            grps.sort()
            for grp in grps:
                lout = "%s\n==Functions in '%s'\b" % (lout,grp)
                lout = "%s\n%s" % (lout,show_list(f[grp]))
        else:
            for nam in args:
                sym = self.tdl.symbolTable.getFunc(nam)
                if sym is None:
                    lout =" cannot find function %s " % nam
                else:
                    lout = "%s:  %s\n" % (sym.name, sym)
                    for i in sym.desc.split('\n'): lout = "%s    %s\n" %(lout,i)
        show_more(lout,writer=self.stdout)

    def show_variable(self,args):
        lout = ''
        if len(args)==0:
            f = self.tdl.symbolTable.listData()
            grps = f.keys()
            grps.sort()
            for grp in grps:
                lout = "%s\n==Variables in '%s'\b" % (lout,grp)
                lout = "%s\n%s" % (lout,show_list(f[grp]))
        else:
            for nam in args:
                sym = self.tdl.symbolTable.getVariable(nam)
                if sym is None:
                    lout =" cannot find variable %s " % nam
                else:
                    lout = " %s:  %s\n" % (sym.name, sym)
                    # for i in sym.desc.split('\n'): lout = "%s    %s\n" %(lout,i)
        show_more(lout,writer=self.stdout)

    def show_groups(self,args=None):
        " print list of groups"
        l = self.tdl.symbolTable.listGroups()
        print "   Default Data     Group = '%s'" % self.tdl.symbolTable.dataGroup
        print "   Default Function Group = '%s'" % self.tdl.symbolTable.funcGroup
        print "   ==Currently defined groups: "
        show_more(show_list(l),writer=self.stdout)

    def show_group(self,grp):
        " list all contents of a group "
        l = self.tdl.symbolTable.listGroups()
        lout = ''
        if grp in l:
            fcns = self.tdl.symbolTable.listFunc()[grp]
            data = self.tdl.symbolTable.listData()[grp]
            if len(fcns)>0:
                lout = "%s\n==Functions in '%s'\b" % (lout,grp)
                lout = "%s\n%s" % (lout,show_list(fcns))
            else:
                lout = "%s\n==No Functions in '%s'\b" % (lout,grp)
            if len(data)>0:
                lout = "%s\n==Variables in '%s'\b" % (lout,grp)
                lout = "%s\n%s" % (lout,show_list(data))
            else:
                lout = "%s\n==No Variables in '%s'\b" % (lout,grp)
            show_more(lout,writer=self.stdout)
        else:
            print " No group %s.  Try 'show groups'" % grp
    def do_shell(self, arg):
        import os
        os.system(arg)

    def do_show(self,argin):
        args = argin.strip().split()
        key = None
        if len(args) > 0:  key = args.pop(0).strip()
        
        if key is None:
            form = "\nDefault Data group = %s\nDefault Func group = %s"
            print form % (self.tdl.symbolTable.dataGroup,self.tdl.symbolTable.funcGroup)
            f = self.tdl.symbolTable.listFunc()
            d = self.tdl.symbolTable.listData()
            for grp in self.tdl.symbolTable.listGroups():
                print "\n===Group: %s" % grp
                for t,name in ((f,'functions'),(d,'variables')):
                    n = 0
                    if t.has_key(grp):  n = len(t[grp])
                    print "   %i %s " % (n,name)
            
        elif key  in ( '-f', 'functions','function'):
            self.show_function(args)

        elif key  in ( '-v', 'variables','variable'):
            self.show_variable(args)

        elif key in ('-g', 'groups'):
            self.show_groups()
            
        elif key in ('group') and len(args)>0:
            self.show_group(args[0].strip())
        else:
            args.insert(0,key)
            self.do_show_symbols(args)

    def do_help(self,argin):
        args = argin.strip().split()
        key = None
        if len(args) > 0:  key = args.pop(0).strip()
        if key is None:
            print self.help('help')
        elif key in ('-t','topics'):
            print "  Additional help is avaiable on the following topics:\n"
            print  show_list(self.help_topics, ncol = 5)
        elif key == 'topic' :
            topic = args[0].strip()
            if topic in self.help_topics:
                show_more(self.help(topic),writer=self.stdout)
            else:
                print "  No help on topic %s. Try 'help topics'" % (topic)
        elif key in self.help_topics:
            show_more(self.help(key),writer=self.stdout)
        else:
            args.insert(0,key)
            self.do_show_symbols(args,msg='help')
            
    def do_show_symbols(self,args,msg='show'):
        for key in args:
            if key.endswith(','): key=key[:-1]
            key.strip()
            if len(key)>0:
                if self.tdl.symbolTable.getFunc(key):
                    self.show_function([key])
                elif self.tdl.symbolTable.getVariable(key):
                    self.show_variable([key])
                elif self.tdl.symbolTable.hasGroup(key):
                    self.show_group(key)
                elif msg == 'show':
                    print " cannot find %s (try 'help show' or 'show groups')" % key
                else:
                    print "  No help on %s.  Try 'help' or 'help topics'" % key

    def tdl_execute(self,s_inp):
        self.default(s_inp)
        
    def default(self,s_inp):
        s = s_inp.strip()
        words = s.split()
        if s in ('quit','exit','EOF'):
            return 1
        elif s.startswith('!'):
            os.system(s)

        else:
            self._status = False
            ret,x,detail = None, None, None
            try:
                if s.startswith('show '):
                    self.do_show(s[5:])
                else:
                    ret = self.tdl.eval(s)
                self._status = True
            except ValueError, detail:
                x = "%s error: %s" % ('syntax',detail) 
            except TypeError, detail:
                x = "%s error: %s" % ('syntax/type',detail) 
            except (ArithmeticError), detail:
                x = "%s error: %s" % ('mathematical',detail) 
            except (LookupError), detail:
                x = "%s error: %s" % ('array lookup',detail) 
            except (EvalException), detail:
                x = "%s error: %s" % ('evaluation',detail) 
            except (ParseException), detail:
                x = "%s error: %s" % ('syntax',detail)            
            except:
                x = 'unknown error'
            if not self._status and x is not None:
                if detail is None: detail = s
                print x
                if self.debug or x == 'unknown error':
                    print '='*60
                    traceback.print_exc()
                    print '='*60

            elif ret is not None:
                print ret
                
    
#####################################################################################               
#####################################################################################
def show_usage():
    print Help.ShellUsage
    sys.exit(1)
    
def main(arg,debug=False):
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hd",["help", "debug"])
    except:
        show_usage()
    for key,val in opts:
        if key in ("-d", "--debug"): debug = True
    t = shell(debug=debug)
    t.cmdloop()

if (__name__ == '__main__'):
    # show_usage()
    main(sys.argv[1:])
    
