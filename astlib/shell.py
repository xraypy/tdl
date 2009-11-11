#!/usr/bin/python2.6
#

import cmd
from os import path, environ, system, getcwd
import sys
import types
import getopt
import traceback
import numpy
import compiler

import inputText
from util import EvalError

try:
    import scipy
    scipy_version = scipy.__version__
except:
    scipy_version = '(not available)'
    

banner = """  Larch %s  M. Newville, T. Trainor (2009)
  using python %s, numpy %s, and scipy %s""" % (compiler.__version__,
                                                '%i.%i.%i' % sys.version_info[:3],
                                                numpy.__version__,
                                                scipy_version)

class shell(cmd.Cmd):
    intro  = "  === Type 'help' to get started ==="
    ps1    = "larch> "
    ps2    = ".....> "
    max_save_lines = 5000
    def __init__(self, completekey='tab', scripts=None, debug=False,
                 stdin=None, stdout=None, intro=None, GUI='TkAgg'):


        print banner
        
        self.debug  = debug
        try:
            import readline
            self.rdline = readline
        except ImportError:
            self.rdline = None

        cmd.Cmd.__init__(self,completekey='tab')

        self.historyfile = path.join(environ.get('HOME',getcwd()),'.larch_history')
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

        
        self.larch  = compiler.Compiler()
        self.input  = inputText.InputText(prompt=self.ps1)
        self.prompt = self.ps1

        
    def __del__(self):
        if (self.rdline):
            self.rdline.set_history_length(1000)
            self.rdline.write_history_file(self.historyfile)

    def emptyline(self):     pass
    def do_shell(self, arg):   system(arg)

    def do_help(self,arg):   self._helpshow(arg, cmd='help')
    def do_show(self,arg):   self._helpshow(arg, cmd='show')

    def larch_execute(self,s_inp):  self.default(s_inp)
        
    def _helpshow(self,arg, cmd='help'):
        if arg.startswith('(') and arg.endswith(')'): arg = arg[1:-1]
        # print 'helpshow ', arg, cmd
        self.default("%s(%s)"% (cmd,repr(arg)))
                     
        
    def default(self,text):
        text = text.strip()
        if text in ('quit','exit','EOF'):
            return 1
        elif text.startswith('!'):
            self.do_shell(text[1:])
        else:
            ret = None
            self.input.put(text,lineno=0)
            while len(self.input) >0:
                block,fname,lineno = self.input.get()
                ret = self.larch.eval(block,fname=fname,lineno=lineno)
                if self.larch.error:
                    i  = self.larch.error[0]                        
                    print "\n".join(i.get_error())
                        
                if ret is not None: print ret

if __name__ == '__main__':
    t = shell(debug=True).cmdloop()
