#!/usr/bin/python2.6
#
import numpy
import compiler

import inputText
from util import EvalError
import cmd
import os
import sys
import types
import getopt
import traceback


class shell(cmd.Cmd):
    banner = """Tiny Data Language %s  M. Newville, T. Trainor (2009)
Using python %s and numpy %s\n"""
    intro  = "   Type 'help' to get started\n"
    ps1    = "tdl> "
    ps2    = "...> "
    max_save_lines = 5000
    def __init__(self, completekey='tab', scripts=None,libs=None, debug=False,
                 stdin=None, stdout=None, intro=None, GUI='TkAgg'):
        self.debug = debug
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

        pyversion = '%i.%i.%i' % sys.version_info[:3]
        print self.banner % (compiler.__version__, pyversion, numpy.__version__)
        
        self.compiler = compiler.Compiler()
        self.eval     = self.compiler.eval
        self.input    = inputText.InputText(prompt=self.ps1,interactive=False)

        self.prompt    = self.ps1
        self._status   = True

                    
    def __del__(self):
        if (self.rdline):
            self.rdline.set_history_length(1000)
            self.rdline.write_history_file(self.historyfile)

    def emptyline(self):     pass
    def do_shell(self, arg):   os.system(arg)

    def do_help(self,arg):   self._helpshow(arg, cmd='help')
    def do_show(self,arg):   self._helpshow(arg, cmd='show')

    def tdl_execute(self,s_inp):  self.default(s_inp)
        
    def _helpshow(self,arg, cmd='help'):
        if arg.startswith("'") and arg.endswith("'"): arg = arg[1:-1]
        if arg.startswith('"') and arg.endswith('"'): arg = arg[1:-1]
        self.default("%s('%s')"% (cmd,arg))
        
    def load_run(self,text,filename=None,lineno=None):
        self.input.put(text,filename=filename, lineno=lineno)
        while len(self.input) >0:
            block,fname,lineo = self.input.get()
            if block is None:
                self.prompt = self.ps2
                return None
            
            ret = self.eval(block)
            self.prompt = self.ps1
            return ret
    
    def default(self,inp):
        s = inp.strip()
        words = s.split()
        if s in ('quit','exit','EOF'):
            return 1
        elif s.startswith('!'):
            os.system(s)
        else:
            self._status = False
            ret,x,detail = None, None, None
            try:
                ret = self.load_run(s)
            except (ValueError,NameError), detail:
                x = "%s error: %s" % ('syntax',detail)
            except TypeError, detail:
                x = "%s error: %s" % ('syntax/type',detail)
            except (ArithmeticError), detail:
                x = "%s error: %s" % ('mathematical',detail)
            except (LookupError), detail:
                x = "%s error: %s" % ('lookup',detail)
            except (EvalError), detail:
                x = "%s error: %s" % ('evaluation',detail)
            except (NameError,AttributeError), detail:
                x = "%s error: %s" % ('evaluation',detail)
            except:
                x = "%s error: %s" % ('syntax',detail)

            if not self._status and x is not None:
                if detail is None: detail = s
                print x
                if self.debug or x == 'unknown error':
                    print '='*60
                    traceback.print_exc()
                    print '='*60
            elif ret is not None:
                try:
                    print ret
                except:
                    print "Error printing ret value"


if (__name__ == '__main__'):
    opts, args = getopt.getopt(sys.argv[1:], "hq", ["help", "quiet"])

    for k,v in opts:
        if k in ("-h", "--help"):   show_usage()
        if k in ("-q", "--quiet"):  quiet = v
        if k in ("-d", "--debug"):  debug = v

    t = shell(debug=debug)
    print args, opts
    if len(args)>0:
        for s in args:
            if s.endswith('.py'): s = s[:-3]
            elif s.endswith('.tdl'): s = s[:-4]
            ret = t.default("import %s" % s)
        ret
    else:
        t.cmdloop()
