#!/usr/bin/env python
############################################################################
#
# -------------
# Modifications
# -------------
# 6-10-06 T2
# - override open, so all open statements use Util.file_open class
#   dangerous???
#
# 4-27-06 T2
# - moved help stuff to Help
#
# 4-8-06 T2
# - small tweak to setting stdin/stdout so tdl_wxGUI will work.
# - small change to do_show_symbols, checks if key is a group.
#
############################################################################

from Num import Num, num_version

from Eval   import Evaluator
from Symbol import SymbolTable
import Util
from Util   import ParseException, EvalException
from Util   import show_list, split_arg_str, show_more
from Util   import PrintExceptErr
from Help   import Help
import version
import cmd
import os
import sys
import types
import getopt
import traceback


class shell(cmd.Cmd):
    banner = """
    Tiny Data Language %s  M. Newville, T. Trainor (2006)
    Using %s\n"""
    intro  = "\n    Type 'help' to get started\n"
    ps1    = "tdl> "
    ps2    = "...> "
    max_save_lines = 500

    def __init__(self, completekey='tab', libs=None, debug=False,
                 stdin=None, stdout=None, intro=None, GUI='TkAgg'):

        cmd.Cmd.__init__(self,completekey='tab')

        self.use_rawinput = True
        if stdin is not None:
            sys.stdin = stdin
            self.use_rawinput = False

        if stdout is not None:  sys.stdout = stdout

        self.stdin = sys.stdin
        self.stdout = sys.stdout

        print self.banner % (version.version,num_version)

        self.tdl       = Evaluator(interactive= True, debug= debug,
                                   input=self.stdin, output = self.stdout,
                                   libs= libs, GUI=GUI)

        self.prompt    = self.ps1
        self.tdl.prompt= self.ps2
        self._status   = True
        #self.debug     = debug

        # override builtin open function
        # is this a bad idea?????
        __builtins__['open'] = Util.file_open(sym=self.tdl)

    def emptyline(self):
        pass

    def do_shell(self, arg):
        import os
        os.system(arg)

    def do_help(self,arg):
        s_inp = "help '%s'" % arg
        self.default(s_inp)

    def do_show(self,arg):
        s_inp = "show '%s'" % arg
        self.default(s_inp)

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
                ret = self.tdl.execute(s)
                self._status = True
            except (ValueError,NameError), detail:
                x = "%s error: %s" % ('syntax',detail)
            except TypeError, detail:
                x = "%s error: %s" % ('syntax/type',detail)
            except (ArithmeticError), detail:
                x = "%s error: %s" % ('mathematical',detail)
            except (LookupError), detail:
                x = "%s error: %s" % ('array lookup',detail)
            except (EvalException,Util.EvalException), detail:
                x = "%s error: %s" % ('evaluation',detail)
            except (ParseException,Util.ParseException), detail:
                x = "%s error: %s" % ('syntax',detail)
            except:
                x = 'unknown error'
            if not self._status and x is not None:
                if detail is None: detail = s
                print x
                if self.tdl.debug or x == 'unknown error':
                    print '='*60
                    traceback.print_exc()
                    print '='*60
            elif ret is not None:
                try:
                    print ret
                except:
                    PrintExceptErr("Error printing ret value")


#####################################################################################
#####################################################################################
#def show_usage():
#    print Help.ShellUsage
#    sys.exit(1)
#
#def main(arg,debug=False):
#    try:
#        opts, args = getopt.getopt(sys.argv[1:], "hd",["help", "debug"])
#    except:
#        show_usage()
#    for key,val in opts:
#        if key in ("-d", "--debug"): debug = True
#    t = shell(debug=debug)
#    t.cmdloop()

if (__name__ == '__main__'):
    # show_usage()
    #main(sys.argv[1:])
    t = shell(debug=debug)
    t.cmdloop()

