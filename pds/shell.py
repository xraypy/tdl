#!/usr/bin/python
##########################################################################
"""
Tom Trainor (fftpt@uaf.edu ), 6-10-2006 
Simple python shell program 

Modifications:
--------------
- modified from tdl-0.2
- major revision 10/08 TPT

"""
##########################################################################
"""
ToDo:

- see wiki

"""
##########################################################################

import os
import sys
import types
import getopt
import time

from   interpretor import Interpretor
import util
from   util import PrintExceptErr, command2expr, show_list
from   util import split_args, show_more, trimstring
from   numimports import _NumShell

##########################################################################

# some global flags
QUIT     = -1   # recieved a quit command
COMPLETE =  0   # normal return --> ps1
CONTINUE =  1   # Continuation --> ps2

# python keywords, make sure a command doesnt clash with these
PYTHON_KEY_WORDS = ['and','as','del','for','is','raise',
                  'assert','elif','from','lambda','return', 
                  'break','else','global','not','None','try',      
                  'class','except','if','or','while'    
                  'continue','exec','import','pass','yield',    
                  'def','finally','in','print']

##########################################################################
HELP_STR =  """\n
***************************** HELP ******************************
===The following help options are available at the command line:
>help              # Help options
>help -u           # General useage notes
>help name         # Provide help on specified object

===To list data and variables, use the show command: 
>show              # List defined data/functions
>show  name        # Detailed list
>help show         # For more show options
*****************************************************************\n
"""

##########################################################################
HELP_USE_STR = """\n
*************************** HELP USE ****************************
PDS is a simple shell wrapper around python - the command syntax
is the same as using the default python interpretor.  The difference
are a few add-ons, some organization and default imports.  For information
about python scripting start at:  http://www.python.org

** Startup / Default Configuration
The main purpose of pds is to make interactive python work easier.
This is in large part achieved by making sure the default namespace
that you are working in at startup has defined all the modules you need
commonly need.
 - more notes on editing startup files

** Commands
The list of commands are dislplayed by typing 'show' at the command line
Commands can be exectuted from the command line as:

  pds>command arg1, arg2

When you exectute a command the command name and arguments get 'caught'
before being sent to the python interpretor.  They are exectuted
in one of two ways, depending on how the command was defined:

  - shell commands are defined in the 'shell' and in general do not send
    the command arguments to the python interpretor.  Examples are
    show, quit, clear, etc..

  - command wrappers are simply command shortcuts to functions, or aliases
    of a complete function/argument call. The command/arg syntax gets
    repacked into a function call and sent to the interpretor.
    Type 'help addcmd' and 'help alias' for more information on how to
    create these.

** Data Organization and Display


** Helpful Commands and Added Builtins


*****************************************************************\n
""" 

##########################################################################
class Shell(_NumShell):
    __doc__  = """
    ****************************************
    * Python Data Shell                    *
    * Type 'help' to get started           *
    ****************************************
    """
    
    max_save_lines = 500

    #############################################################
    def __init__(self,args=[],stdin=None,stdout=None,
                 completekey='tab',intro=None,debug=False,GUI='TkAgg'):
        # Prompt
        self.ps1 = "pds>"
        self.ps2 = "...>"
        
        # Set stdin/out
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

        # Startup strings
        self.intro = self.__doc__
        if intro:
            self.intro = self.intro + intro
        if self.intro:
            self.stdout.write(str(self.intro)+"\n")

        # Interpretor
        self.interp = Interpretor()

        # Data and vars
        self.debug         = debug
        self.GUI           = GUI
        self.queue         = []
        self.commands      = {}
        self.file_error_break = True

        # Builtin commands
        self.pds_commands  = {'quit':self.do_quit,
                              'exit':self.do_quit,
                              'EOF':self.do_quit,
                              'show':self.do_show,
                              'help':self.do_help,
                              'load':self.do_load,
                              'debug':self.do_debug,
                              'execfile':self.do_execfile,
                              'save':self.do_save,
                              'restore':self.do_restore,
                              'clear':self.do_clear,
                              'addcmd':self.do_addcmd,
                              'alias':self.do_addalias}

        # run startup 
        self.startup(args=args)

        # check
        for cmd in self.pds_commands.keys():
            if cmd in PYTHON_KEY_WORDS:
                print "Warning: command '%s' is a python keyword" % cmd

        # setup readline if have it. -- add later
        self.completekey = completekey

    #############################################################
    def startup(self,args=[]):
        """
        set up builtins and exec startup arguments
        """
        
        # Builtins
        from pds.builtins import __pdsbuiltins__
        self.interp.symbol_table.data["__pds__"] = __pdsbuiltins__
        
        # add module functions to __builtins__
        startup = []
        startup.append("__builtins__.update({'group':__pds__.group})")
        startup.append("__builtins__.update({'ls':__pds__._ls})")
        startup.append("__builtins__.update({'cd':__pds__._cd})")
        startup.append("__builtins__.update({'pwd':__pds__._cwd})")
        startup.append("__builtins__.update({'mod_import':__pds__.mod_import})")
        startup.append("__builtins__.update({'path':__pds__._path})")

        # functions we want to use with command syntax
        self.do_addcmd('ls',"__pds__.ls")
        self.do_addcmd('pwd',"__pds__.pwd")
        self.do_addcmd('cd',"__pds__.cd")
        self.do_addcmd('more',"__pds__.more")
        self.do_addcmd('path',"__pds__.path")
        self.do_addcmd('mod_import',"__pds__.mod_import")

        # if numeric/scientific stuff is included
        # load it. Note leave this here in case there
        # are cmd line args that depend on numerics...
        if hasattr(self,'num_setup'):
            tmp = self.num_setup()
            startup = startup + tmp
        
        # any other args add to queue
        if type(args) == types.ListType:
            for a in args:
                s = a.strip()
                startup.append(s)
        self.queue = startup

        return

    ############################################################################
    def loop(self, intro=None):
        """
        Repeatedly issue a prompt, accept input, parse an initial prefix
        off the received input, and dispatch to action methods, passing them
        the remainder of the line as argument.  Modified from cmd.Cmd.cmdloop.

        Note for the python interpretor to work it is important that blank
        lines are passed, and that lines are not stripped entireley
        """
        stop = None
        while stop != QUIT:
            if self.queue:
                line = self.queue.pop(0)
            else:
                if stop == COMPLETE:
                    prompt = self.ps1
                else:
                    prompt = self.ps2
                if self.use_rawinput:
                    try:
                        line = raw_input(prompt)
                    except EOFError:
                        line = 'EOF'
                else:
                    #line = raw_input(self.prompt)
                    self.stdout.write(prompt)
                    self.stdout.flush()
                    line = self.stdin.readline()
                    #if not len(line):
                    if line == None:
                        line = 'EOF'
                    else:
                        line = line[:-1] # chop \n
            # exec the line
            stop = self.exec_line(line)
    
    #############################################################
    def exec_line(self, line):
        """
        This method first checks for blank line.
        If not blank then
            If a shell command: pass to shell
        Inspect the first token of command line and checks:
            If special quit keywords return quit: (quit, exit, EOF)
            If in list of command keywords: repack as function and exec
        Otherwise exec via interpretor

        Note for the python interpretor to work it is important that blank
        lines are passed, and that lines are not stripped entireley
        
        """
        s = str(line).rstrip()
        try:
            # parse
            if len(s) > 0:
                if s.startswith('!'):
                    os.system(s[1:])
                    return COMPLETE
                
                words = s.split()
                cmd = words[0].strip()
                if len(words) > 1:
                    idx = s.find(cmd)+len(cmd)
                    arg = s[idx:].strip()
                else:
                    arg = ''
                if self.debug: print '**cmd: %s  %s' % (cmd,arg)
            else:
                cmd = None
            
            # if its a cmd execute it
            if cmd:
                if cmd in self.pds_commands.keys():
                    return self.exec_cmd(cmd,arg)
                elif cmd in self.commands.keys():
                    return self.exec_cmd(cmd,arg)
            
            # otherwise pass to the interpretor
            return self.interp.execute(s)
            
        except:
            err_str = "Error executing line:\n%s" % s
            if self.debug:
                PrintExceptErr(err_str,print_trace=True)
            else:
                PrintExceptErr(err_str,print_trace=False)
            return COMPLETE

    #############################################################
    def exec_cmd(self,cmd,arg):
        """
        exectute commands
        """
        ret = None
        if cmd in self.pds_commands.keys():
            fun = self.pds_commands[cmd]
            ret = fun(arg)
        elif cmd in self.commands.keys():
            cmd_fun = self.commands[cmd]
            s = cmd_fun + ' ' + arg
            s = command2expr(s,symtable=self.interp.symbol_table)
            if self.debug: print '**command repacked: ', s
            ret = self.interp.execute(s)
        else:
            print "Error: uknown command: %s" % cmd
        if ret == None:
            return COMPLETE
        else:
            return ret

    #############################################################

    #############################################################
    ## Builtin commands
    #############################################################

    #############################################################
    def do_quit(self,arg):
        """
        Quit the shell 
        """
        try:
            self.close_pylab()
        except:
            pass
        return QUIT
    
    #############################################################
    def do_debug(self,arg):
        """
        Toggle debug flag 
        """
        if self.debug == True:
            print "Debug off"
            self.debug = False
        else:
            print "Debug on"
            self.debug = True

    #############################################################
    def do_load(self,fname):
        """
        Load file of pds commands for execution
        Note that file execution will halt if there is an
        error in the file.  This default behaviour can
        be changed by setting the flag
        file_error_break = False
        """
        if os.path.exists(fname) and os.path.isfile(fname):
            f = open(fname)
            lines = f.readlines()
            f.close()
        else:
            print 'File error: cannot find file to load for %s ' % fname
            return
        
        for line in lines:
            try:
                stop = self.exec_line(line)
                if stop == QUIT: return QUIT 
            except:
                if self.file_error_break == True:
                    err_str =  "Exception in file %s" % fname
                    util.PrintExceptErr(err_str)
                    return COMPLETE
        if self.debug: print 'load done'
        return COMPLETE

    #############################################################
    def do_execfile(self,fname):
        """
        Execute a file of python code
        """
        return self.interp.execute_file(fname)

    #############################################################
    def do_addcmd(self,*arg):
        """
        Allows a function to be called with command syntax.
          >>addcmd mycmd, myfun
        results in
          >>mycmd arg1, arg2, key=xx
        being repackaged to
          >>myfun(arg1,arg2,key=xx)
        which is then sent to the interpretor
        """
        if len(arg) == 1:
            arg = arg[0]
            words = arg.split(',')
            if len(words) != 2:
                print "Error parsing cmd/func names"
            cmd_name = words[0].strip()
            func_name = words[1].strip()
        elif len(arg) == 2:
            cmd_name  = arg[0].strip()
            func_name = arg[1].strip()
        else:
            print "Wrong number of arguments"
            return
        cmd_name  = trimstring(cmd_name)
        func_name = trimstring(func_name)

        if len(cmd_name) == 0 or len(func_name) == 0:
            print "Error parsing cmd/func names"
            return
        
        # do a check
        if cmd_name in PYTHON_KEY_WORDS:
            print "Error: command '%s' is a python keyword" % cmd_name
            return
        
        if cmd_name in self.pds_commands.keys():
            print "Error: command '%s' is a builtin keyword" % cmd_name
            return

        if cmd_name in self.commands.keys():
            print "Warning: overwriting command '%s' " % cmd_name
            self.commands[cmd_name] = func_name
        else:
            self.commands.update({cmd_name:func_name})
        return

    #############################################################
    def do_addalias(self,*arg):
        """
        Create a command shortcut for a function.
          >>alias myalias, 'myfun(args)'
        results in
          >>myalias
        being repackaged to
          >>myfun(args)
        which is then sent to the interpretor

        Note alias's are added to the command list
        (see 'show' and 'help addcmd')

        example
        pds>alias "code", "cd('~/code')"
        pds>code
        pds>pwd
        /home/bob/code
        """
        if len(arg) == 1:
            arg = arg[0]
            words = arg.split(',')
            if len(words) != 2:
                print "Error parsing cmd/func names"
            cmd_name  = words[0].strip()
            func_name = words[1].strip()
        elif len(arg) == 2:
            cmd_name  = arg[0].strip()
            func_name = arg[1].strip()
        else:
            print "Wrong number of arguments"
            return
        cmd_name  = trimstring(cmd_name)
        func_name = trimstring(func_name)

        if len(cmd_name) == 0 or len(func_name) == 0:
            print "Error parsing cmd/func names"
            return

        # do a check
        if cmd_name in PYTHON_KEY_WORDS:
            print "Error: command '%s' is a python keyword" % cmd_name
            return

        s = "__%s__ = lambda : %s" % (cmd_name,func_name)
        ret = self.exec_line(s)
        s = "addcmd '%s', '__%s__'" % (cmd_name,cmd_name)
        ret = self.exec_line(s)

    ##############################################################

    ##############################################################
    ## save / restore state and clear variables
    ##############################################################
    
    ##############################################################
    def do_save(self,args):
        """
        Save program state to a file
        Note this may fail since not all
        objects can be pickled... needs improvement!
        Call as:
        >save            # save all to default fname
        >save fname      # save all to fname
        >save data fname # save data to file
        """
        from util import pickle_1 as pickle
        #from util import pickle_2 as pickle

        # parse input, get filename
        args = split_args(args)
        dname = None
        if len(args) == 0:
            t = time.strftime("%Y_%m_%d_%H%M", time.localtime())
            fname = 'save_%s.sav' % t
        else:
            if len(args) == 1:
                fname = args[0]
            elif len(args) == 2:
                dname = args[0]
                fname = args[1]
            else:
                return
            
        # get data dictionary
        data = self.interp.symbol_table.get_data_dict(name=dname)
        if data == None: return

        # pickle it    
        pickle(data,fname)
        
    ##############################################################
    def do_restore(self,fname):
        """
        Restore state from a file
        """
        from util import unpickle_1 as unpickle
        #from util import unpickle_2 as unpickle

        pdata = unpickle(fname)
        if pdata == None:
            print "No data"
            return
        if len(pdata) > 0:
            self.interp.symbol_table.put_data_dict(pdata)
            
    ##############################################################
    def do_clear(self,arg):
        """
        Clear all 'data' from the workspace
        """
        args  = split_args(arg)
        if len(args) == 0:
            dname = [None]
        else:
            dname = args
        for d in dname:
            data = self.interp.symbol_table.get_data_dict(name=d)
            for key in data.keys():
                del self.interp.symbol_table.data[key]

    ##############################################################

    #############################################################
    ##  Help and Show
    #############################################################
                
    ##############################################################
    def do_help(self,arg):
        """
        \rhelp: display help.
        \rrun '>help' to display options
        """
        if arg == None or len(arg) == 0:
            #print self.help_str
            print HELP_STR
            return
        
        (opts, args) = getopt.getopt(arg.split(), "u",["use"])
        for key,val in opts:
            if key in ("-u", "--use"):
                #print self.help_use_str
                print HELP_USE_STR
        for a in args:
            if a in self.pds_commands.keys():
                #tmp = getattr(self,a)
                tmp = self.pds_commands[a]
                help(tmp)
            elif a in self.commands.keys():
                cmd_fun = self.commands[a]
                s = "help(%s)" % cmd_fun
                self.interp.execute(s)
            else:
                s = "help(%s)" % a
                self.interp.execute(s)
        return

    #############################################################
    def do_show(self,arg):
        """
        \rshow:  list functions and variables.
        \r Options:
        \r >show -ath group
        \r  -a = show all (include symbols with '_' in list)
        \r  -b = show builtins (defined in __builtins__ dictionary).  Note this option
        \r       trumps tunnel and args (ie -t and other args are ignored).  
        \r  -t = tunnel, display attributes within a class instance or module.
        \r       Note modules are only displayed with one level of tunneling
        \r       (ie -t option is ignored if a module is passed as an argument).
        \r       However, this will tunnel on all class instances.
        """
        tunnel=False
        skip=True
        (opts, args) = getopt.getopt(arg.split(), "abt",["all","builtins","tunnel"])
        for key,val in opts:
            if key in ("-a", "--all"):
                skip=False
            if key in ("-b", "--builtin"):
                self._show_builtins(_skip=skip)
                return
            if key in ("-t", "--tunnel"):
                tunnel=True
        if len(args)>0:
            #arg = ' '.join(args)
            for a in args:
                self._show(symbol=a,tunnel=tunnel,_skip=skip)
        else:
            self._show(tunnel=tunnel,_skip=skip)

    #############################################################
    def _show_builtins(self,_skip=True):
        d = self.interp.symbol_table.list_builtins(_skip=_skip)
        print "\n***** Builtins ******"
        self._print_show(d)

    #############################################################
    def _show(self,symbol=None,tunnel=False,_skip=True):
        if symbol:
            ty = self.interp.symbol_table.sym_type(symbol)
            if ty == None:
                print "Symbol '%s' not found" % symbol
                return
            elif ty in ('oth'):
                print "Symbol '%s' of uknown type" % symbol
                print type(self.interp.symbol_table.get_symbol(symbol))
                return
        else: 
            ty = None
            print "\n==== Commands ===="
            tmp = self.commands.keys()
            tmp = tmp + self.pds_commands.keys()
            tmp.sort()
            print show_list(tmp,textwidth=82)

        if ty == 'var':
            tmp = self.interp.symbol_table.get_symbol(symbol)
            print repr(tmp)
            return
        elif ty == 'fnc':
            tmp = self.interp.symbol_table.get_symbol(symbol)
            print help(tmp)
            return
        else:
            d = self.interp.symbol_table.list_symbols(symbol=symbol,
                                                    tunnel=tunnel,
                                                    _skip=_skip)
            self._print_show(d)

    #############################################################
    def _print_show(self,d):
        # Print stuff
        if d == None: return
        if len(d['mod'])> 0:
            print "\n==== Modules ===="
            print show_list(d['mod'],textwidth=82)
        if len(d['fnc'])> 0:
            print "\n==== Functions / Classes ===="
            print show_list(d['fnc'],textwidth=82)
        if len(d['ins'])> 0:
            print "\n==== Instances ===="
            print show_list(d['ins'],textwidth=82)
        if len(d['var'])> 0:
            print "\n==== Variables ===="
            print show_list(d['var'],textwidth=82)
        if len(d['oth'])> 0:
            print "\n==== Other Python Objects ===="
            print show_list(d['oth'],textwidth=82)
        print "\n"
        return

    ##############################################################

##################################################################################### 
#####################################################################################
def show_usage():
    print 'Startup options for pds:'
    print '-h:  help'
    print '-w:  use wx GUI'
    print '-d:  debug on'
    print '-v:  verbose'
    print '-x:  non-interactive execution'
    print 'Examples (combine the various options):'
    print '>pds -v                  #verbose startup'
    print '>pds -vd                 #verbose startup and debug on'
    print '>pds -w                  #use wx GUI'
    print '>pds -v  fname1 fname2   #verbose and execute files'
    print '>pds fname1 fname2       #defaults and execute files'
    sys.exit(1)
    
#####################################################################################
def main(arg):
    ##############################################################
    # test for command line switches
    ##############################################################
    use_wxGui = False
    debug     = False
    verbose   = False
    do_shell  = True
    files     = []

    opts, args = getopt.getopt(arg, "wdvhxl:")
    for o, a in opts:
        if o == '-w':
            use_wxGui = True
        elif o == '-d':
            debug = True
        elif o == '-v':
            verbose = True
        elif o == '-x':
            do_shell = False
        elif o == '-h':
            show_usage()
        
    #################################################################
    # Import and set default paths:
    #   pds path is the path of the pds directory
    #   lib_path is the parent directory of pds
    #   mods_path is either lib_path/modules or lib_path
    # Additional paths should be set in the startup files
    #################################################################
    tmp = globals()
    __path__  = os.path.abspath(tmp['__file__'])
    pds_path  = os.path.dirname(__path__)
    lib_path  = os.path.split(pds_path)[0]
    mods_path = os.path.join(lib_path,"modules")
    if not os.path.exists(mods_path):
        mods_path = lib_path

    # Set Python Path
    if verbose:
        print ' == pds paths:'
        print '    pds_path   = %s  ' % pds_path
        print '    lib_path   = %s  ' % lib_path
        if os.path.exists(mods_path):
            print '    mods_path  = %s  ' % mods_path
    util.set_path('.',recurse=False,verbose=verbose)
    util.set_path(pds_path,recurse=False,verbose=verbose)
    util.set_path(lib_path,recurse=False,verbose=verbose)
    util.set_path(mods_path,recurse=True,verbose=verbose)

    ##############################################################
    # startup  files:
    #   first site-wide ==> pds_path/startup.pds
    #   then from user's home dir ==> ~/.pds
    # files contains tuples of (file, Warn_if_Not_Exist_Flag)
    ##############################################################
    files = [(os.path.join(pds_path,'startup.pds'),False)]
    user_home = os.path.expanduser('~')
    files.append((os.path.join(user_home,'.pds'),False))

    # the rest of the arg line are added files
    if args is not None:
        for f in args: files.append((f,True))

    if verbose:
        print " == Startup files:"
        for f in files: print f

    # create a dictionary of default system variables
    sys_vars = {}
    sys_vars['__pds__.home'] = user_home
    args = []
    for var in sys_vars.keys():
        #shell.interp.setVariable(var,sys_vars[var])
        s = "%s='%s'" % (var,sys_vars[var])
        args.append(s)

    ##########################################################
    # run it
    ##########################################################
    if use_wxGui == False:
        s = Shell(args=args,debug=debug,GUI='TkAgg')
        for f,warn in files:
            if os.path.exists(f) and os.path.isfile(f):
                s.do_load(f)
            elif warn:
                print "\n  ***Cannot find file: %s" % f
        if do_shell:
            s.loop()
    elif use_wxGui == True:
        # pds gets started from within the wxGui
        # looks like dir gets reset when call application
        work_dir = os.getcwd()
        from wxgui import wxShell
        wxShell.intro     = None
        wxShell.debug     = debug
        wxShell.files     = files
        wxShell.args      = args
        rsrc_path         = os.path.join(mods_path,'wxgui')
        wxShell.rsrc_path = rsrc_path
        rsrc = os.path.join(rsrc_path,'wxShell.rsrc.py')
        gui  = wxShell.model.Application(wxShell.wxShell, aFileName=rsrc)
        os.chdir(work_dir)
        gui.MainLoop()

#####################################################################################
#####################################################################################
if (__name__ == '__main__'):
    # show_usage()
    main(sys.argv[1:])
