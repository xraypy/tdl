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
############################################################################

from Num import Num, num_version

from Eval   import Evaluator
from Symbol import SymbolTable
import Util
from Util   import show_list, split_arg_str, show_more
from Util   import EvalError, ParseError, SymbolError, ConstantError
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

    def __init__(self, completekey='tab', scripts=None,libs=None, debug=False,
                 stdin=None, stdout=None, intro=None, GUI='TkAgg'):
        try:
            import readline
            self.rdline = readline
        except ImportError:
            self.rdline = None

        cmd.Cmd.__init__(self,completekey='tab')

        self.historyfile = os.path.join(os.environ.get('HOME',os.getcwd()),'.tdl_history')
        if self.rdline is not None:
            try:
                self.rdline.read_history_file(self.historyfile)
            except IOError:
                pass
            
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
        # __builtins__['open'] = Util.file_open(sym=self.tdl)
        if scripts is not None:
            for i in scripts:
                self.default('load("%s")' % i)
                
    def __del__(self):
        
        if (self.rdline):
            self.rdline.set_history_length(1000)
            self.rdline.write_history_file(self.historyfile)

    def emptyline(self):
        pass

    def do_shell(self, arg):
        import os
        os.system(arg)

    def _helpshow(self,arg, cmd='help'):
        if arg.startswith("'") and arg.endswith("'"): arg = arg[1:-1]
        if arg.startswith('"') and arg.endswith('"'): arg = arg[1:-1]
        self.default("%s('%s')"% (cmd,arg))
        
    def do_help(self,arg):   self._helpshow(arg, cmd='help')
    def do_show(self,arg):   self._helpshow(arg, cmd='show')


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
                # print " tdl execute : " , s
                ret = self.tdl.execute(s)
                self._status = True
            except (ParseError,Util.ParseError), detail:
                x = "%s error: %s" % ('syntax',detail)
            except (ValueError,NameError), detail:
                x = "%s error: %s" % ('syntax',detail)
            except TypeError, detail:
                x = "%s error: %s" % ('syntax/type',detail)
            except (ArithmeticError), detail:
                x = "%s error: %s" % ('mathematical',detail)
            except (LookupError), detail:
                x = "%s error: %s" % ('array lookup',detail)
            except (EvalError,Util.EvalError), detail:
                x = "%s error: %s" % ('evaluation',detail)
            except (SymbolError,AttributeError,Util.ParseError), detail:
                x = "%s error: %s" % ('evaluation',detail)
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
                    self.tdl.ShowError("Error printing ret value")


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
    debug=False

    t = shell(debug=debug,     scripts = sys.argv[1:])

    t.cmdloop()

