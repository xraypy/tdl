"""
Simple python shell program 

Authors / Modifications:
------------------------
* Tom Trainor (tptrainor@alaska.edu ), 6-10-2006 
* modified from tdl-0.2: Matt Newville (newvile@cars.uchicago.edu)
* major revision 10/08 TPT

Notes:
------
To start the pds shell program use the runpds.py script:
(this script makes sure that the root directory of pds is
on the sys.path so imports will work).  
>>python runpds.py

Or make sure that the pds directory is on your path and:
>>python
>>from pds import shell
>>shell.main()

More information can be found at:
http://cars9.uchicago.edu/ifeffit/tdl/Docs/Pds

"""
##########################################################################

import os
import sys
import types
import getopt
import time
try:
    import readline
except:
    pass

from  .interpretor import Interpretor
from  .shellutil   import set_path, show_more, show_list
from  .shellutil   import PrintExceptErr, command2expr
from  .shellutil   import split_args, trimstring, split_cmd_line
from  .numshell    import _NumShell

##########################################################################

# some global flags
QUIT       = -1  # recieved a quit command
SUCCESS    =  0  # all ok
COMPLETE   =  0  # normal return --> ps1
CONTINUE   =  1  # continuation --> ps2
EXECERROR  =  2  # there was an error executing a line
TEXT_WIDTH = 62

# python keywords, make sure a command doesnt clash with these
PYTHON_KEY_WORDS = ['and','as','del','for','is','raise',
                    'assert','elif','from','lambda','return', 
                    'break','else','global','not','None','try',
                    'class','except','if','or','while',    
                    'continue','exec','import','pass','yield',
                    'def','finally','in','print']
PYTHON_KEY_WORDS.sort()

##########################################################################
HELP_STR =  """\n
***************************** HELP ******************************
===The following help options are available at the command line:
>>help              # Help options
>>help -u           # General useage notes
>>help name         # Provide help on specified object

===To list data and variables, use the show command: 
>>show              # List defined data/functions
>>show name         # Detailed list
>>help show         # For more show options
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
that you are working in at startup has defined all the modules you 
commonly need.

This can be customized via the startup files:
   INSTALL_PATH/pds/startup.pds
     - and -
   HOMEPATH/.pds

You can determine the values for these paths using the following:
pds>show __pds__.__path__   # INSTALL_PATH
pds>show __home__           # HOMEPATH

** Syntax
Command line syntax is the same as at the typical python interpretor prompt
with some minor differences.  We have defined a list of 'commands' that can
be called with a command syntax (see below).  You can also call the shell
using the following syntax for example: '>>!dir'.  Finally modules and functions
etc can be displayed in an organized format using the '>>show' command.  

** Commands
The list of commands (and other loaded objects) are displayed by typing
'show' at the command line. Commands can be exectuted from the
command line as:

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

More information on pds and the overall project can be found at:
http://cars9.uchicago.edu/ifeffit/tdl/Docs/Pds

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
        """
        Init
        """
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
        self.error_break   = True
        self.commands      = {}

        # Builtin commands
        self.pds_commands  = {'quit':self.do_quit,
                              'exit':self.do_quit,
                              #'EOF':self.do_quit,
                              'show':self.do_show,
                              'help':self.do_help,
                              'load':self.do_load,
                              'debug':self.do_debug,
                              'fexec':self.do_execfile,
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
        from builtins import __pdsbuiltins__
        self.interp.symbol_table.data["__pds__"] = __pdsbuiltins__
        
        # add module functions to __builtins__
        startup = []
        startup.append("__builtins__.update({'group':__pds__.group})")
        startup.append("__builtins__.update({'ls':__pds__._ls})")
        startup.append("__builtins__.update({'cd':__pds__._cd})")
        startup.append("__builtins__.update({'pwd':__pds__._cwd})")
        startup.append("__builtins__.update({'rimport':__pds__.rimport})")
        startup.append("__builtins__.update({'path':__pds__._path})")
        startup.append("__builtins__.update({'source':__pds__.source})")
        startup.append("__builtins__.update({'info':__pds__.info})")

        # functions we want to use with command syntax
        self.do_addcmd('ls',"__pds__.ls")
        self.do_addcmd('pwd',"__pds__.pwd")
        self.do_addcmd('cd',"__pds__.cd")
        self.do_addcmd('more',"__pds__.more")
        self.do_addcmd('path',"__pds__.path")
        self.do_addcmd('rimport',"__pds__.rimport")
        self.do_addcmd('source',"__pds__.source")
        self.do_addcmd('info',"__pds__.info")

        # if numeric/scientific stuff is included
        # load it. Note leave this here in case there
        # are cmd line args below that depend on numerics...
        if hasattr(self,'num_setup'):
            tmp = self.num_setup()
            startup = startup + tmp
        
        # any other args add to queue
        if type(args) == types.ListType:
            for a in args:
                s = a.strip()
                startup.append(s)
        self.queue = startup

        # execute the queue here before moving on
        # otherwise startup stuff is not available
        # to files executed before we start the loop
        # also allow it to run through even if there is
        # an error...
        self.error_break = False
        self.exec_queue()
        self.error_break = True

        return

    ############################################################################
    def loop(self,):
        """
        The prompt loop
        
        Repeatedly issue a prompt, accept input, parse an initial prefix
        off the received input, and dispatch to action methods, passing them
        the remainder of the line as argument.  Modified from cmd.Cmd.cmdloop.

        Note for the python interpretor to work it is important that blank
        lines are passed, and that lines are not stripped entirely (ie leave
        white space on left)
        """
        stop = COMPLETE
        while stop != QUIT:
            if stop == CONTINUE:
                prompt = self.ps2
            else:
                prompt = self.ps1
            if self.use_rawinput:
                try:
                    line = raw_input(prompt)
                except EOFError:
                    line = 'quit'
            else:
                self.stdout.write(prompt)
                self.stdout.flush()
                line = self.stdin.readline()
                if line == None:
                    line = 'quit'
                else:
                    #line = line[:-1]    # chop \n
                    line = line.rstrip() # chop \n
            # split on ';' and exec the line
            lines = split_cmd_line(line)
            for s in lines:
                stop = self.exec_line(s)
                if stop == QUIT: break
        return SUCCESS
    
    #############################################################
    def exec_queue(self):
        """
        Execute lines in self.queue
        """
        if len(self.queue) == 0:
            return  SUCCESS
        for line in self.queue:
            # split on ';' and exec the line
            lines = split_cmd_line(line)
            for s in lines:
                stop = self.exec_line(s)
                if stop == QUIT: 
                    self.queue = []
                    return QUIT
                elif (stop == EXECERROR) and (self.error_break == True):
                    self.queue = []
                    print "Error executing shell.queue"
                    return EXECERROR
        return SUCCESS
    
    #############################################################
    def exec_line(self, line):
        """
        Exec a single line
        
        This method first inspect the first token of the command line
            If its a '!' send it to the shell
            If its a 'command' call _exec_cmd()
            Otherwise pass to the interpretor

        The return value from this method should always be one of:
            QUIT, SUCCESS, COMPLETE, CONTINUE, EXECERROR
        Exceptions from the interpretor are handled here.
        
        Note for the python interpretor to work it is important that blank
        lines are passed in and that lines are not stripped entirely
        (ie leave left white space). We also let the interpretor handle
        comments (#).
        """
        s = str(line).rstrip()
        try:
            if len(s) > 0:
                #if self.debug: print '**line:%s' % (s)
                # see if its shell cmd
                if s.startswith('!'):
                    os.system(s[1:])
                    return COMPLETE
                # look for commands
                words = s.split()
                cmd = words[0].strip()
                if len(words) > 1:
                    idx = s.find(cmd)+len(cmd)
                    arg = s[idx:].strip()
                else:
                    arg = ''
            else:
                cmd = None
            # if its a cmd execute it in _exec_cmd
            if cmd:
                if cmd in self.pds_commands.keys():
                    return self._exec_cmd(cmd,arg)
                elif cmd in self.commands.keys():
                    return self._exec_cmd(cmd,arg)
            # otherwise pass to the interpretor
            return self.interp.execute(s)
        except:
            err_str = "Error executing line:\n'%s'" % s
            if self.debug:
                PrintExceptErr(err_str,print_trace=True)
            else:
                #PrintExceptErr(err_str,print_trace=False)
                PrintExceptErr(err_str,print_trace=True)
            # clear the exception and interpretors buffer
            # sys.exc_clear()
            self.interp.console.resetbuffer()
            return EXECERROR
        
    #############################################################
    def _exec_cmd(self,cmd,arg):
        """
        execute commands
        
        Note only call this from exec_line for proper
        error / exception handling! 
        """
        ret = None
        # local commands
        if cmd in self.pds_commands.keys():
            fun = self.pds_commands[cmd]
            ret = fun(arg)
        # other commands
        elif cmd in self.commands.keys():
            cmd_fun = self.commands[cmd]
            s = cmd_fun + ' ' + arg
            s = command2expr(s,symtable=self.interp.symbol_table)
            if self.debug: print '**command repacked: ', s
            ret = self.interp.execute(s)
        else:
            print "Error: uknown command: %s" % cmd
        #if ret == None:
        if ret not in (QUIT, SUCCESS, COMPLETE, CONTINUE, EXECERROR):
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
            self.close_pyplot()
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
        return SUCCESS

    #############################################################
    def do_load(self,fname):
        """
        Load file of pds commands for execution

        Note that file execution will halt if there is an
        error in the file.  This default behaviour can
        be changed by setting the flag self.error_break = False

        If the file has a '.sav' extension we first try to read it
        using restore (ie assume its a pickle)
        """
        if fname[-4:] == '.sav':
            if self._restore(fname) == 1: return COMPLETE
        
        if os.path.exists(fname) and os.path.isfile(fname):
            f = open(fname)
            lines = f.readlines()
            f.close()
        else:
            print 'File error: cannot find file to load for %s ' % fname
            return COMPLETE
        if len(lines) == 0: return COMPLETE

        # execute the lines through the queue
        self.queue = []
        self.queue.extend(lines)
        ret = self.exec_queue()
        if ret == EXECERROR:
            print "Error in file %s" % fname
        if self.debug: print 'load done'
        return COMPLETE

    #############################################################
    def do_execfile(self,fname):
        """
        Execute a file of python code

        Note that the file will be executed as if it was
        imported into the '__main__' namespace
        """
        if not os.path.exists(fname):
            files = []
            for p in sys.path:
                f = os.path.join(p,fname)
                if os.path.exists(f):
                    files.append(f)
            if len(files) == 0:
                print "File '%s' not found" % fname
                return SUCCESS
            elif len(files) == 1:
                fname = files[0]
            else:
                print "Warning multiple files found on path"
                print "Please specify the full path in fname or change the"
                print "working directory to the directory with the correct file"
                for f in files: print "    %s" % f
                return SUCCESS
        return self.interp.execute_file(fname)

    #############################################################
    def do_addcmd(self,*arg):
        """
        Add a command interface to a function
        
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
            return SUCCESS
        cmd_name  = trimstring(cmd_name)
        func_name = trimstring(func_name)

        if len(cmd_name) == 0 or len(func_name) == 0:
            print "Error parsing cmd/func names"
            return SUCCESS
        
        # do a check
        if cmd_name in PYTHON_KEY_WORDS:
            print "Error: command '%s' is a python keyword" % cmd_name
            return SUCCESS
        
        if cmd_name in self.pds_commands.keys():
            print "Error: command '%s' is a builtin keyword" % cmd_name
            return SUCCESS

        if cmd_name in self.commands.keys():
            print "Warning: overwriting command '%s' " % cmd_name
            self.commands[cmd_name] = func_name
        else:
            self.commands.update({cmd_name:func_name})
        return SUCCESS

    #############################################################
    def do_addalias(self,*arg):
        """
        Create a command shortcut for a function.

        The alias is just a way of making shortcuts for repetitive taks
          >>alias myalias, 'myfun(args)'
        results in
          >>myalias
        being repackaged to
          >>myfun(args)
        which is then sent to the interpretor

        Note alias's are added to the command list
        (see 'show' and 'help addcmd')

        Example:
        --------
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
            return SUCCESS
        cmd_name  = trimstring(cmd_name)
        func_name = trimstring(func_name)

        if len(cmd_name) == 0 or len(func_name) == 0:
            print "Error parsing cmd/func names"
            return SUCCESS

        # do a check
        if cmd_name in PYTHON_KEY_WORDS:
            print "Error: command '%s' is a python keyword" % cmd_name
            return SUCCESS

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

        Examples:
        ---------
        >>save            # save all to default fname
        >>save fname      # save all to fname
        >>save data fname # save data to file
        """
        from pds.shellutil import pickle_1 as pickle
        #from pds.shellutil import pickle_2 as pickle

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
        if data == None: return SUCCESS

        # pickle it    
        pickle(data,fname)
        
    ##############################################################
    def do_restore(self,fname):
        """
        Restore state from a file that was created with save
        """
        self._restore(fname)
        return SUCCESS 

    def _restore(self,fname):
        from pds.shellutil import unpickle_1 as unpickle
        #from pds.shellutil import unpickle_2 as unpickle

        pdata = unpickle(fname)
        if pdata == None:
            print "No data"
            return 0
        if len(pdata) > 0:
            self.interp.symbol_table.put_data_dict(pdata)
            return 1
        return 0
    
    ##############################################################
    def do_clear(self,arg):
        """
        Clear all 'data' from the workspace

        This clears everything that looks like a variable, ie
        simple data types and class instances
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
        return SUCCESS

    ##############################################################

    #############################################################
    ##  Help and Show
    #############################################################
    
    ##############################################################
    def do_help(self,arg):
        """
        Display help.

        Example:
        --------
        >>help     # displays help options
        """
        if arg == None or len(arg) == 0:
            #print self.help_str
            print HELP_STR
            return SUCCESS
        
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
                ret = self.exec_line(s)
            else:
                s = "help(%s)" % a
                ret = self.exec_line(s)
        return SUCCESS

    #############################################################
    def do_show(self,arg):
        """
        List functions and variables.

        Options:  >>show -ath group
        --------
        -a = show all (include symbols with names that have a leading '_' )
        -b = show builtins (defined in __builtins__ ).  Note this option
             trumps tunnel (ie -t is ignored).  
        -t = tunnel, display attributes within a class instance or module.
             Note modules are only displayed with one level of tunneling
             (ie -t option is ignored if a module is passed as an argument).
             However, this will tunnel on all class instances.
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
        return SUCCESS

    #############################################################
    def _show_builtins(self,_skip=True):
        """
        """
        print "\n***** Builtins ******"
        # show python keywords
        print "\n==== Python Key Words ===="
        print show_list(PYTHON_KEY_WORDS,textwidth=TEXT_WIDTH)
        # get builtins
        d = self.interp.symbol_table.list_builtins(_skip=_skip)
        self._print_show(d)
        return SUCCESS

    #############################################################
    def _show(self,symbol=None,tunnel=False,_skip=True):
        """
        """
        if symbol:
            ty = self.interp.symbol_table.sym_type(symbol)
            if ty == None:
                print "Symbol '%s' not found" % symbol
                return SUCCESS
            elif ty in ('oth'):
                print "Symbol '%s' of uknown type" % symbol
                print type(self.interp.symbol_table.get_symbol(symbol))
                return SUCCESS
        else: 
            ty = None
            print "\n==== Commands ===="
            tmp = self.commands.keys()
            tmp = tmp + self.pds_commands.keys()
            tmp.sort()
            print show_list(tmp,textwidth=TEXT_WIDTH)

        if ty == 'var':
            tmp = self.interp.symbol_table.get_symbol(symbol)
            print repr(tmp)
            return SUCCESS
        elif ty == 'fnc':
            tmp = self.interp.symbol_table.get_symbol(symbol)
            print help(tmp)
            return SUCCESS
        else:
            d = self.interp.symbol_table.list_symbols(symbol=symbol,
                                                      tunnel=tunnel,
                                                      _skip=_skip)
            self._print_show(d)
        return SUCCESS

    #############################################################
    def _print_show(self,d):
        """
        Print stuff
        """
        if d == None: return SUCCESS
        if len(d['mod'])> 0:
            print "\n==== Modules ===="
            print show_list(d['mod'],textwidth=TEXT_WIDTH)
        if len(d['fnc'])> 0:
            print "\n==== Functions / Classes ===="
            print show_list(d['fnc'],textwidth=TEXT_WIDTH)
        if len(d['ins'])> 0:
            print "\n==== Instances ===="
            print show_list(d['ins'],textwidth=TEXT_WIDTH)
        if len(d['var'])> 0:
            print "\n==== Variables ===="
            print show_list(d['var'],textwidth=TEXT_WIDTH)
        if len(d['oth'])> 0:
            print "\n==== Other Python Objects ===="
            print show_list(d['oth'],textwidth=TEXT_WIDTH)
        print "\n"
        return SUCCESS

##################################################################################### 
#####################################################################################
def show_usage():
    """
    Print use options
    """
    print 'Startup options for pds:'
    print '-h:  help'
    print '-w:  use wx GUI'
    print '-d:  debug on'
    print '-v:  verbose'
    print '-x:  non-interactive execution'
    print 'Examples (combining the various options):'
    print '>pds -v                  #verbose startup'
    print '>pds -vd                 #verbose startup and debug on'
    print '>pds -w                  #use wx GUI'
    print '>pds -v  fname1 fname2   #verbose and execute files'
    print '>pds fname1 fname2       #defaults and execute files'
    sys.exit(1)
    
#####################################################################################
def main(arg=''):
    """
    Startup the shell program
    """
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
    # Default paths. Assume the following layout:
    #    pds_path  = root1/root2/pds
    #    startup   = root1/root2/pds/startup.pds
    #
    # Here we try to make sure that the PYTHONPATH includes:
    #    ".", "root1", "root1/root2", and "root1/root2/pds"
    #
    # Additional paths should be set in the startup files
    #################################################################
    #tmp = globals()
    #__path__  = os.path.abspath(tmp['__file__'])
    pds_path  = os.path.split(os.path.abspath(__file__))[0]
    root2     = os.path.split(pds_path)[0]
    root1     = os.path.split(root2)[0]

    # Set Python Path
    if verbose:
        print ' == pds paths:'
        print '    root1  = %s  ' % root1
        print '    root2  = %s  ' % root2
        print '    pds_path   = %s  ' % pds_path
    set_path('.',recurse=False,verbose=verbose)
    set_path(root2,recurse=False,verbose=verbose)
    set_path(root1,recurse=False,verbose=verbose)
    set_path(pds_path,recurse=False,verbose=verbose)

    ##############################################################
    # startup  files:
    #   first site-wide ==> pds_path/startup.pds
    #   then from user's home dir ==> ~/.pds
    # files contains tuples of (file, Warn_If_Not_Exist_Flag)
    ##############################################################
    files = [(os.path.join(pds_path,'startup.pds'),False)]
    user_home = os.path.expanduser('~')
    files.append((os.path.join(user_home,'.pds'),False))

    # the rest of the arg line are added files
    if args is not None:
        for f in args: files.append((f,True))

    if verbose:
        print " == Startup files:"
        for f in files: print '    ', f

    # create a dictionary of default system variables
    sys_vars = {}
    if sys.platform == 'win32':
        pds_path  = pds_path.replace('\\','/')
        root_path = root2.replace('\\','/')
        user_home = user_home.replace('\\','/')
    sys_vars['__pds__.__path__']     = pds_path
    sys_vars['__pds__.__rootpath__'] = root_path
    sys_vars['__home__']             = user_home
    args = []
    for var in sys_vars.keys():
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
        from pds.pcgui import wxShell
        wxShell.intro     = None
        wxShell.debug     = debug
        wxShell.files     = files
        wxShell.args      = args
        rsrc_path         = os.path.join(pds_path,'pcgui')
        wxShell.rsrc_path = rsrc_path
        rsrc = os.path.join(rsrc_path,'wxShell.rsrc.py')
        gui  = wxShell.model.Application(wxShell.wxShell, aFileName=rsrc)
        os.chdir(work_dir)
        gui.MainLoop()

        
def rungui():
    """
    Startup the shell program
    """
    ##############################################################
    # test for command line switches
    ##############################################################
    use_wxGui = True
    debug     = False
    verbose   = False
    do_shell  = True
    files     = []
    args  = None
    
    #################################################################
    # Default paths. Assume the following layout:
    #    pds_path  = root1/root2/pds
    #    startup   = root1/root2/pds/startup.pds
    #
    # Here we try to make sure that the PYTHONPATH includes:
    #    ".", "root1", "root1/root2", and "root1/root2/pds"
    #
    # Additional paths should be set in the startup files
    #################################################################
    #tmp = globals()
    #__path__  = os.path.abspath(tmp['__file__'])
    pds_path  = os.path.split(os.path.abspath(__file__))[0]
    root2     = os.path.split(pds_path)[0]
    root1     = os.path.split(root2)[0]

    # Set Python Path
    if verbose:
        print ' == pds paths:'
        print '    root1  = %s  ' % root1
        print '    root2  = %s  ' % root2
        print '    pds_path   = %s  ' % pds_path
    set_path('.',recurse=False,verbose=verbose)
    set_path(root2,recurse=False,verbose=verbose)
    set_path(root1,recurse=False,verbose=verbose)
    set_path(pds_path,recurse=False,verbose=verbose)

    ##############################################################
    # startup  files:
    #   first site-wide ==> pds_path/startup.pds
    #   then from user's home dir ==> ~/.pds
    # files contains tuples of (file, Warn_If_Not_Exist_Flag)
    ##############################################################
    files = [(os.path.join(pds_path,'startup.pds'),False)]
    user_home = os.path.expanduser('~')
    files.append((os.path.join(user_home,'.pds'),False))

    # the rest of the arg line are added files
    if args is not None:
        for f in args: files.append((f,True))

    if verbose:
        print " == Startup files:"
        for f in files: print '    ', f

    # create a dictionary of default system variables
    sys_vars = {}
    if sys.platform == 'win32':
        pds_path  = pds_path.replace('\\','/')
        root_path = root2.replace('\\','/')
        user_home = user_home.replace('\\','/')
    sys_vars['__pds__.__path__']     = pds_path
    sys_vars['__pds__.__rootpath__'] = root_path
    sys_vars['__home__']             = user_home
    args = []
    for var in sys_vars.keys():
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
        from pds.pcgui import wxShell
        wxShell.intro     = None
        wxShell.debug     = debug
        wxShell.files     = files
        wxShell.args      = args
        rsrc_path         = os.path.join(pds_path,'pcgui')
        wxShell.rsrc_path = rsrc_path
        rsrc = os.path.join(rsrc_path,'wxShell.rsrc.py')
        gui  = wxShell.model.Application(wxShell.wxShell, aFileName=rsrc)
        os.chdir(work_dir)
        gui.MainLoop()

#####################################################################################
