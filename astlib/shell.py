#!/usr/bin/python2.6
#

import cmd
import os
import sys
import types
import getopt
import traceback
import numpy

import compiler
import inputText
from util import EvalError

class shell(cmd.Cmd):
    banner = """Tiny Data Language %s  M. Newville, T. Trainor (2009)
    with python %s and numpy %s"""
    intro  = "  //=== Type 'help' to get started ==="
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
        
        self.compiler = compiler.Compiler(load_builtins=True)
        self.input    = inputText.InputText(prompt=self.ps1)

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
        
    def default(self,text):
        text = text.strip()
        if text in ('quit','exit','EOF'):
            return 1
        elif text.startswith('!'):
            os.system(text[1:] )
        else:
            self._status = False
            ret,x,detail = None, None, None
            try:
                self.input.put(text)
                while len(self.input) >0:
                    block,fname,lineno = self.input.get()
                    ret = self.compiler.eval(block)

            except (ValueError,NameError), detail:
                x = "Shell %s error: %s" % ('syntax',detail)
            except TypeError, detail:
                x = "Shell %s error: %s" % ('syntax/type',detail)
            except (ArithmeticError), detail:
                x = "Shell %s error: %s" % ('mathematical',detail)
            except (LookupError), detail:
                x = "Shell %s error: %s" % ('lookup',detail)
            except (EvalError), detail:
                x = "Shell %s error: %s" % ('evaluation',detail)
            except (NameError,AttributeError), detail:
                x = "Shell %s error: %s" % ('evaluation',detail)
            except:
                x = "Shell %s error: %s" % ('syntax',detail)

            if not self._status and x is not None:
                if detail is None: detail = s
                print x
                print 'FNAME ', fname, lineno
                print sys.exc_info()
                ast = self.compiler.compile(block)
                print self.compiler.dump(ast,annotate_fields=True,
                                         include_attributes=True)
                if self.debug or x == 'unknown error':
                    print '=*'*40
                    # traceback.print_exc()
                    print '=*'*40
            elif ret is not None:
                try:
                    print ret
                except:
                    print "Error printing ret value"

if (__name__ == '__main__'):
    t = shell(debug=True).cmdloop()
