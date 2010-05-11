#!/usr/bin/env python2.6
#
from __future__ import print_function
import cmd
import os
import sys
import numpy

import interpreter
import inputText

HISTFILE = '.larch_history'
BANNER = """  Larch %s  M. Newville, T. Trainor (2009)
  using python %s, numpy %s"""


class shell(cmd.Cmd):
    ps1    = "larch> "
    ps2    = ".....> "
    max_save_lines = 5000
    def __init__(self, completekey='tab', scripts=None, debug=False,
                 stdin=None, stdout=None, quiet=False,userbanner=None):

        if not quiet:
            print(BANNER % (interpreter.__version__,
                            '%i.%i.%i' % sys.version_info[:3],
                            numpy.__version__))

            if userbanner is not None:
                print(userbanner)
            
        self.debug  = debug
        try:
            import readline
            self.rdline = readline
        except ImportError:
            self.rdline = None
        cmd.Cmd.__init__(self,completekey='tab')
        homedir = os.environ.get('HOME', os.getcwd())

        self.historyfile = os.path.join(homedir, HISTFILE)
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
        
        self.larch  = interpreter.Interpreter()
        self.input  = inputText.InputText(prompt=self.ps1)
        self.prompt = self.ps1
        
    def __del__(self):
        if (self.rdline):
            self.rdline.set_history_length(1000)
            self.rdline.write_history_file(self.historyfile)

    def emptyline(self):
        pass

    def parseline(self, line):
        """Parse the line into a command name and a string containing
        the arguments.  Returns a tuple containing (command, args, line).
        'command' and 'args' may be None if the line couldn't be parsed.
        """
        line = line.strip()
        if not line:
            return None, None, line
        elif line[0] == '?':
            line = 'help ' + line[1:]
        elif line[0] == '!':
            if hasattr(self, 'do_shell'):
                line = 'shell ' + line[1:]
            else:
                return None, None, line
        return '', '', line

    def do_shell(self, arg):
        os.system(arg)

    def larch_execute(self,s_inp):
        self.default(s_inp)

    def loadfile(self,filename):
        fh = open(filename,'r')
        for i,line in enumerate(fh.readlines()):
            self.input.put(line, filename=filename,lineno=i)

    def default(self,text):
        text = text.strip()
        if text in ('quit','exit','EOF'):
            return 1

        if text.startswith('help'):
            arg = text[4:]
            if arg.startswith('(') and arg.endswith(')'): arg = arg[1:-1]
            if arg.startswith("'") and arg.endswith("'"): arg = arg[1:-1]
            if arg.startswith('"') and arg.endswith('"'): arg = arg[1:-1]
            text  = "help(%s)"% (repr(arg))
        if text.startswith('!'):
            self.do_shell(text[1:])
        else:
            ret = None
            self.input.put(text,lineno=0)
            
            self.prompt = self.ps2
            while len(self.input) >0:
                block,fname,lineno = self.input.get()
                ret = self.larch.eval(block,fname=fname,lineno=lineno)
                if callable(ret) and not isinstance(ret,type):
                    try:
                        if 1 == len(block.split()):
                            ret = ret()
                    except:
                        pass
                if self.larch.error:
                    err = self.larch.error.pop(0)
                    fname, lineno = err.fname, err.lineno
                    print("%s: %s" % err.get_error())

                    for err in self.larch.error:
                        if ((err.fname != fname or err.lineno != lineno)
                            and err.lineno >= 0):
                            print('%s' % (err.get_error()[1]))

                elif ret is not None:
                    print("%s" % ret)
                self.prompt = self.ps1
            
if __name__ == '__main__':
    t = shell(debug=True).cmdloop()
