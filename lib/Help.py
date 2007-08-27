# Help class for on-line documentation /help in tdl
# M. Newville Univ of Chicago, T. Trainor, Univ of Alaska Fairbanks (2005,2006)
#
# --------------
# Modifications
# --------------
#  28-Apr-06 T2 - moved show/help from shell into Help class
#
#  20-Mar-06 MN - changed to have a Help class that contains a dictionary of
#                 help topics
#                
#
##########################################################################

ShellUsage = """
  Tdl command-line shell help:
   
"""

HelpHelp =  """
  ###########################################################################
   = Help options =
    help                  # Help options (this list)
    help topics  (-t)     # list topics for additional help
    help <name>           # Detailed help on specific topic, function, or variable
    (see also: 'help show')
  ###########################################################################
"""

HelpShow =  """
  ###########################################################################
   = show: show list groups and variables =
    show              # summary of all top-level groups
    show -l           # short list of all top-level groups
    show mygroup      # list variables and sub-groups of mygroup
    show myvariable   # summary for variable myvariable
    (see also: 'help' and 'help topics')
  ###########################################################################
"""

Overview = """
This is TDL (Tiny Data Language, Matt Newville and Tom Trainor 2006).

Tdl is a data processing language written in python and designed to be 
  - easy to use for novices.
  - complete enough for intermediate to advanced data processing.
  - data-centric, so that arrays of data are easy to manage and use.
  - easily extendable with python.
 
Tdl uses Numeric Python and SciPy for all its numerical functionality, so
processing arrays of numerical data is easy and fast, and a wide range
of 'advanced numeric capabilities' such as FFTs and non-linear least squares
minimization is (potentially) available for Tdl.

TDL Features:

==Simple Syntax:
  Variable and function names are simple, undecorated strings.

==Data Types:
  Variabls can contain any of the following data types:  ints, floats,
    complex, strings, lists, dictionaries (hashes), or NumPy arrays.

  Multi-dimensional data can be accessed either element-by-element or with
  more complicated "slices":
        tdl> x   = 1             # int
        tdl> y   = 7 + 2j        # complex number:  real + imag'j'
        tdl> arr = arange(10)    # simple array 
        tdl> print arr 
            [0 1 2 3 4 5 6 7 8 9 10]
        tdl> print arr[3]
            3
        tdl> print arr[2:8]
            [2 3 4 5 6 7]
        tdl> str = 'here is a string'
        tdl> print str[0:6]
            here i

==Full Featured, extensible:  many useful functions for numerical data
  processing are built in to Tdl.  More importantly, it is easy to add more
  functions, either from Python or as Tdl procedures (see below). 

==Namespaces:  Each variable and function has two parts to its name:
  a Group Name and Variable Name.  The fully qualified name includes a '.'
  and looks like group.var.

  Unqualified names (ie, a name without a 'group.' prefix) will use a default
  order to determine which Group to associate with a name:  There is always
  a "Default Data Group" and a "Default Function Group" each of which can be
  changed during a TDL session.  At startup, there are two pre-defined groups:
  "_builtin", which contains many built in functions and constants (sin and pi,
  for example), and "_main".  At startup, "_main" is empty, and is set as the
  "Default Data Group" and "Default Function Group".
  
  When creating or assigning a variable without an explicit group name, the default
  Data Group is used.  When creating or assigning a procedure without an explicit
  group name, the default Function Group is used.
  
  When using (say, on the right-hand side of a statement) a variable or function with
  an unqualified group name, the first matching name in the list
     (Default Data Group, Default Function Group, "_main", "_builtin") 

  is used.  The Default Data Group can be changed with
        tdl> datagroup('mydata')

  and the Default Function Group can be changed with 
        tdl> funcgroup('myfuncs')

  Assigning a fully qualified variable or function to a group that previously did
  not exist will create that group:
        tdl> x = 1
        tdl> group1.x = 3.3
        tdl> print x
            1
        tdl> datagroup('group1')
        tdl> print x
            3.3
        tdl> print _main.x, group1.x
            1 3.3
        tdl> show groups
            Default Data Group = group1
            Default Function Group = _main
            All Groups:
                _main  _builtin  group1
  
==Functions and Commands.  Functions and commands have similiar name structure as
  data (group.func,group.cmd). 
  The default group name for functions and commands is _main (functions and commands with
  this group name do not need full name qualification).
  Eg. you can also assign variables variable values with setvar() function:
        tdl> setvar('v', 2.3, group='group1')
        tdl> print group1.v

  The potential advantage here is that 'v' and 'group1' are string types,
  not literal variable names, so they can be computed or passed in a
  procedure. 

==Clean syntax for programming:
  The syntax is "Python Inspired", with some notable exceptions.  Most
  importantly, indentation level is not significant, and blocks end with
  with explicit 'end' statements.  A typical for-block will look like this:

        for i in range(10):
            print i
        endfor     

  There are also while blocks and if-elif-else blocks:
        n = 0
        while n<10: 
           print ' n = ',n
           if n==3: print 'here is 3!!'
               n = n + 1
           endwhile
     
        if n > 10:
            print 'No'
        elif n > 5 and n < 8:
            print 'Maybe'
        else:
            print 'Yep'
        endif

  A design goal is that well-formed tdl code should be very easy to
  translate into valid python (and vice versa).

==User-defined functions (aka procedures):
  User defined functions can be written in tdl:
        def  myfunc(arg1, option='test'):
             'documentation string'
              print 'this is my funcition ', arg1
              print 'optional param = ', option
              if type(option) != 'string':     
                  print 'option must be a string!!'
                  return False
              endif
              value = sqrt(arg1)
              x.tmp = value
              return value > 10.
        enddef

  which could be called as 
        tdl> ret = myfunc(3., option = 'xx')
        tdl> print ret, x.tmp
          False 1.73205080757

==dofile() function:
  you can run a file of tdl commands with the dofile() function:
        tdl> dofile('myfile.tdl')

==eval() function:
  you can construct a tdl expression on the fly from strings and execute it:

        tdl> eval("%s.value = %f" %  ['group', 10.2])
        tdl> print group.value 
          10.2

==read_ascii() function:
  you can read in ASCII column data files very easily.

        tdl> read_ascii('my.dat', group='f')
          ['f', 'x', ''y']
      
  this read 2 columns (labeled 'x' and 'y') from the column file and
  created the array variables f.x and f.y.  Also created was f.titles
  to hold the titles in the data file (non-numeric content at the top of
  the file) and f.column_labels

==On-line help:
  well, this is in progress....

==Easy to add your python functions, including getting access to all  
  the 'data groups' inside your python function.

==Differences between tdl and Python (for python users):
  -  tdl has many builtins and assumes Numerical data. 
  -  tdl has no tuples. Lists are used in their place.
  -  indentation does not matter. blocks are ended with 'end***' 
     (if / endif , for /endfor, def/enddef)
  -  when in doubt, assignment makes a copy, and does not give a reference
  -  data types are NOT objects, and so have no methods.  You must use
     procedural approaches 
       x = arange(100)
       reshape(x,[10,10])   instead of x.shape = (10,10)
"""

HelpStrings = """
  Working with strings:

  strings are sequence of text characters. To specify a string, you enclose it in quotes,
  either single or double quotes:
      tdl>  x = 'this is a string'
      tdl>  x = "A string with a ' in it "

  There are a few special characters that can be included in strings by 'escaping' them with a
  backlash character "\\".  The most important of these are '\\n' to write a newline character,
  and '\\t' to write a tab character.  Other escaped characters are:
      \\n   newline
      \\r   carriage return
      \\f   form feed
      \\v   vertical tab
      \\b   backspace
      \\t   tab
      \\a   bell
      \\"   literal "
      \\'   literal '
      \\\   the backslach character itself.
      
  You can escape quote characters ('\\"' or '\\'') so that they do not mark the end of a
  string sequence: 
      tdl> x = 'a string\\'s string'
      tdl> y = "\\"a string with a double quote\\", he said."

 
  ==Multi-line Strings with Triple Quotes
   
  Tdl allows multi-line strings (that is strings that span lines).  One simple way to do this is
  to include a newline character ('\\n') character in the string, but this is not always sufficient.
  Another way is to use triple quotes (three quotes in a row: either ''' or  \"\"\") to enclose the
  string.  With this approach, the newlines of the enclosed string are preserved:

     tdl> long_string = '''Here is a mult-line string
     ...> line #2
     ...> line #3 of the long string '''

     tdl> print long_string
     Here is a mult-line string
     line #2
     line #3 of the long string

  As with single quote strings, you have the choice of using single or double quotes, and can
  use escaped character and escaped quotes.

  ==String Formatting

  It is often desirable to create strings from program data. To do this, you *format a string*.
  This is done by putting format codes in a string and supplying values to fill in.  Format
  codes use the '%' character to mark the formatting to do.

  Following normal conventions, a formatted string for a floating point number might look like this: 

      tdl> print "x = %8.4f" % sqrt(12)
      x =   3.4641
      
  The "%8.4f" tells the formatting to format a floating point number (the 'f') with 8 total numbers
  and 4 numbers after the decimal point.   A plain '%' is used to separate the format string and the
  value(s) to format.  Other format codes:

      ......
  
  ==Using a dictionary for string formatting

    Borrowing from Python, tdl also allow you to format a string using a dictionary instead of a
    simple list.  This can be a big advantage when formatting many values or formatting data from
    complicated data structures.  In this method, you explicitly set the dictionary key to use by
    naming them between the '%' and the format code:
    
      tdl> data  = {'name':'File 1', 'date':'Monday, April 3, 2006', 'x': 12.4}
      tdl> print " %(name)s , x = %(x)f" % data
      File 1 , x = 12.4
     
"""

HelpDicts = """
  Working with dictionaries:

"""

HelpArrays = """
  Working with arrays:

"""

HelpLists = """
  Working with lists:
  
"""

HelpControl = """
   Help on Programming Control structures (conditionals, loops, etc)

   Program Control structures are important for controlling what parts of
   a program are run.  If you're familiar with other programming languages,
   the concepts and syntax here should be fairly straightforward.

   The basic concept here is to define blocks of code (that is multiple lines
   of code) to be run as a group.  This allows the block of code to be run
   only under some conditions, or to be run repeatedly with different values
   for some variables in the block.
   
   The syntax of tdl uses a colon ':' and a small number of keywords to define
   the blocks of code:
       if x==0:
          print 'cannot divide by zero'
       else:
          x = 2 /x
       endif
   

 = If statements:

   if statements run a block of code only if some condition is met.  Since this
   is so common, there are a few variations allowed in tdl.  First is the
   one line version:

       if x<0:  x = -x

   which will set x to -x (that is, this will make sure x>=0).  In general, the
   if statement looks like:
       if condition:

   where the ending ':' is important (to indicate where the condition end).
   The "condition" is a statement that is evaluated as a boolean value (either
   True or False -- all numeric values are True except 0.0).  The boolean
   operations 'and', 'or' , and 'not' to construct complex conditions:

       if x>0 and x<10:

       if (x>0 and (b>1 or c>1)):


   A warning: in many programming languages multiple conditions can be relied upon
   to only evaluate until the truth of the statement can be determined.  This is NOT
   the case in tdl: the entire statement may be evaluated.

   The second version uses the same 'if condition:', but on a line by itself, followed
   by  multiple statements.  This version ends with 'endif' to indicate how far the 
   block of code extends:
       if x < 0  :
           x = -x
           print ' x set to ',x 
       endif

   Next, there's the 'if-else-endif' version which runs either one block or the
   other, and looks like this (note the 
       if x<0:
           x = -x
           print ' x set to ',x 
       else:
           print ' x was positive to begin with'
       endif

   The final and most complete variation uses 'elif' (pronounced "else if") to allow
   multiple conditions to be tested:
   
       if x<0:
           x = -x
           print ' x set to ',x 
       elif x>10:
           x = 10
           print ' x was too big, set to ', x
       else:
           print ' x was OK to begin with'
       endif

   Multiple 'elif' blocks can be given, though all 'elif' statements must be
   after the 'if' statements, and an 'else' statement (if present) must be last.
   And the entire construct must end with an 'endif'

   In such constructs, the first block with a condition that is True, and no other
   blocks will be run.

 = For loops:

   for loops repeat a block of code, usually incrementing some value used in
   the loop.  Usually, a for loop runs a predictable number of times (see the
   section below on Break and Continue for when it does not!). A for loop
   looks like:
    
       for i in [1,2,3]:
           print i, i*i
       endfor

   The basic syntax for the first line is 'for <variable> in <list>:'  The list
   can be either a real list or a numerical array -- a very common approach is to
   use the range() function to generate a list of numbers:
       for i in range(10):
           print i, 10-i, i*(10-i)
           if i < 2: print ' i < 2!! '
       endfor

   The block is run repeatedly with the loop variable i set to each value in the
   list.

   There is also a 'one-line' version of the for loop:
  
       for i in [1,2,3]: call_some_function(i)
   

 = While loops:

   while loops are similar to for loops, repeatedly running a block of code
   as long as some condition is true:

      x = 1
      while x<4:
         print x
         x = x + 1
      endwhile

   prints
     1.0
     2.0
     3.0

   Beware that a while loop makes it easy to make an "infinite loop" (for example, if
   the value of x had not been increased in the block).
   
 = Break and Continue in For and While loops

   In both for loops and while loops, the block can be exited with a 'break' statement.
   
   
 = Try / Except blocks:

   Sometimes error happen, and it's desirable to run some block of code

       x = 0
       a = 2.
       try:
           y = a/ x
       except:
           print "OK,  that didn't work so well!"
       endtry
      
"""
#  one quote mark: ' to fix emacs colorization
#########################################################
import types, sys
from Symbol import isGroup, isSymbol, symTypes
from Util   import show_list, split_arg_str, show_more, SymbolError

class Help:
    """ basic help mechanism for Tdl """
    #charlen = 16
    charlen = 24
    def __init__(self,tdl=None,output=None):
        self.tdl = tdl
        self.output = output
        if self.output is None: self.output = sys.stdout

        self.buff = []
        h = self.topics = {}
        h['help']  = HelpHelp
        h['show']  = HelpShow
        h['shell_usage'] = ShellUsage
        h['overview'] = Overview
        h['strings'] = HelpStrings
        h['dicts'] =   HelpDicts
        h['arrays'] =  HelpArrays
        h['lists'] =   HelpLists
        h['control'] = HelpControl
        self.help_topics = self.list_topics()

    def list_topics(self):
        x = self.topics.keys()
        x.sort()
        return x

    def add_topic(self,name,text):
        if type(name) != types.StringType:
            name = str(name)
        self.topics[name] = text

    def get_help(self,name):
        if name in self.topics.keys():
            return self.topics[name]
        return "No help available for %s" % (name)

    ##############################################
    def bprint(self,x):
        self.buff.append("%s\n" % x)
        return

    def __showbuff(self):
        show_more(self.buff)
        self.buff = []
        
    def show(self,arg=None, extended=False):
        if arg is None or arg=='':
            return self.show_topgroups()
        try:
            sym = self.tdl.symbolTable.getSymbol(arg)
            if  isSymbol(sym):
                self.show_symbols([sym], extended=extended)
            elif isGroup(sym):
                self.show_symbols(sym.values(),extended=extended)
            #elif type(arg)==types.FunctionType:
            #    self.bprint(arg.__doc__)
            #    self.__showbuff()
            elif type(arg)==types.StringType:
                args = arg.replace('=',' ').strip().split()
                self.show_symbols(args)
        except:
            print "Show failed"

    def help(self,argin):
        args = argin.strip().split()
        key = None
        if len(args) > 0:  key = args.pop(0).strip()
        self.__showbuff()
        self.help_topics = self.list_topics()
        if key is None:
            self.bprint( self.topics['help'])
        elif key in ('-t','topics'):
            self.bprint("\n==Additional help is avaiable on the following topics:")
            self.bprint(show_list(self.help_topics, ncol = 5))
            self.bprint('')
        elif key in self.help_topics:
            self.bprint(self.get_help(key))
        else:
            args.insert(0,key)
            for a in args:
                self.show(a, extended=True)
        self.__showbuff()
        
    def show_symbols(self,args,indent=1,title=None,groupname=None,followgroups=False,extended=False):
        " list all contents of a group "

        # print "==Show Symbols " , args, indent, title, groupname
        vtab = '  '*indent
        ttab = '  '*(indent-1)
        if title is not None: self.bprint("%s==%s==" % (ttab,title))
        if groupname is None: groupname=''
        args = self.sort_symbols(args)
        for sym in args:
            nam  = sym.name
            if nam.startswith("%s."%groupname): nam = nam[len(groupname)+1:]
            nam  = nam + ' '*(self.charlen - len(nam))
            #print ' >> ', nam,sym
            #self.bprint("%s%s: %s" % (vtab,nam, sym.getinfo(extended=True)))
            inf = sym.getinfo(extended=extended)
            txt = "%s%s: %s" % (vtab,nam,inf)
            self.bprint(txt)
            if isGroup(sym) and followgroups:
                #print ':: ', sym.keys(), sym.values()
                self.show_symbols(sym.values(),
                                  title=sym.getinfo(extended=extended),
                                  groupname=sym.name, 
                                  indent=indent+1,
                                  followgroups=followgroups)
        self.__showbuff()
        
    def show_topgroups(self):
        stable   = self.tdl.symbolTable
        locgroup = stable.LocalGroup
        modgroup = stable.ModuleGroup

        self.bprint("Default Data group = %s,  Function group = %s\n" % (locgroup, modgroup))

        groups = stable.data.values()
        #groups.sort()
        groups = self.sort_symbols(groups)
        self.show_symbols(groups, extended=True)

        self.bprint(" ")
        gnames = [locgroup]
        if modgroup not in gnames: gnames.append(modgroup)
        for g in gnames:
            grp = stable.getSymbol(g)
            if len(grp.keys())>0:
                self.show_symbols(grp.values(),
                                  indent=2,
                                  title="%s: %s" % (g,grp.getinfo(extended=False)),
                                  groupname=grp.name)
        print " "

    def show_topgroups_fmt(self):
        stable   = self.tdl.symbolTable
        locgroup = stable.LocalGroup
        modgroup = stable.ModuleGroup
        sgroups  = stable.searchGroups

        print "=== Default Data Group = %s,  Function Group = %s" % (locgroup, modgroup)
        print "    Search Groups: ", sgroups
        print " "
        
        print "=== Symbols Defined in Default Data Group:"
        ncol = 8
        txt = {'g':[],'f':[],'v':[]}
        try:
            for sym in stable.getSymbol(locgroup).values():
                name = sym.name
                if name.startswith("%s."%locgroup):
                    name = name[len(locgroup)+1:]
                if isGroup(sym):
                    txt['g'].append(name)
                elif sym.type in (symTypes.defpro,symTypes.pyfunc):
                    txt['f'].append(name)
                else:
                    txt['v'].append(name)
            txt['g'].sort()
            txt['f'].sort()
            txt['v'].sort()
            if len(txt['g']) > 0:
                t = "  Subgroups:  "
                cnt = 0
                for g in txt['g']:
                    t = "%s %s,  " % (t,g)
                    cnt = cnt + 1
                    if cnt > ncol:
                        t = "%s \n              " % t
                        cnt = 0
                print t
            if len(txt['f']) > 0:
                t =  "  Functions:  "
                cnt = 0
                for f in txt['f']:
                    t = "%s %s,  " % (t,f)
                    cnt = cnt + 1
                    if cnt > ncol:
                        t = "%s \n              " % t
                        cnt = 0
                print t
            if len(txt['v']) > 0:
                t = "  Variables:  "
                cnt = 0
                for v in txt['v']:
                    t = "%s %s,  " % (t,v)
                    cnt = cnt + 1
                    if cnt > ncol:
                        t = "%s \n              " % t
                        cnt = 0
                print t
            print " "
        except:
            print "          Cannot read default group"

        print "=== Top Level Groups"        
        ncol = 4
        groups = stable.data.values()
        groups = self.sort_symbols(groups)
        txt = {'g':[],'f':[],'v':[],'s':[]}
        ngrps = 0
        for grp in groups:
            if isGroup(grp): 
                txt['g'].append("* Group: %s" % (grp.name))
                (nvar,nfun,ngrp) = grp.stats()
                txt['f'].append(" %i functions" % (nfun))
                txt['v'].append(" %i variables" % (nvar))
                txt['s'].append(" %i sub groups" % (ngrp))
                ngrps = ngrps +1

        j = 0; k = 0; 
        while True:
            t = ''
            for k in range(ncol):
                if j+k < ngrps:
                    t = "%s %-20s\t" % (t,txt['g'][j+k])
            print t
            t= ''
            for k in range(ncol):
                if j+k < ngrps:
                    t = "%s %-20s\t" % (t,txt['f'][j+k])
            print t
            t = ''
            for k in range(ncol):
                if j+k < ngrps:
                    t = "%s %-20s\t" % (t,txt['v'][j+k])
            print t
            t = ''
            for k in range(ncol):
                if j+k < ngrps:
                    t = "%s %-20s\t" % (t,txt['s'][j+k])
            print t
            ##
            print "\n"
            j = j + k + 1
            if j >= ngrps: break
        return

    def sort_symbols(self, syms):
        """ Given a list of symbols, return the list sorted by name"""
        #syms.sort()
        try:
            # force sort
            syms2 = []
            syms2.append(syms.pop(0))
            while len(syms) > 0:
                j = 0
                while j < len(syms2):
                    if syms2[j].name > syms[0].name:
                        break
                    j = j+1
                syms2.insert(j,syms.pop(0))
            return syms2
        except:
            return syms

###########################################################################################
# 
#     def _group(g,indent=0):
#         if not isGroup(g): return None
#         tab = '   '*indent
#         for k,v in g.items():
#             print "%s %s %s" %(tab, k,  v.getinfo())
#             if isGroup(v): show_group(v,indent=indent+1)
# 
#         
#     def _groupshow(self,group, groupname='',indent=1):
#         vtab = '  '*(indent+1)
#         print 'This is _groupshow ', group.name, indent, group, group.keys()
#         
#         self.__showbuff()
#         gname = group.name
#         if gname.startswith(groupname) and gname!=groupname:
#             gname = gname[len(groupname):]
#         for sym in group.values():
#             nam  = sym.name
#             if nam.startswith(gname):  nam = nam[len(gname)+1:]                    
#             if isGroup(sym):
#                 self.bprint("%s%s: group %s" %(vtab,nam,sym.getinfo()))
#                 print ':: ', sym, sym.name
#                 print self.groupstats.keys()
#                 self._groupshow(sym,indent=indent+1, groupname=sym.name)
#             elif isSymbol(sym):
#                 nam  = nam + ' '*(16-len(nam))
#                 self.bprint( "%s%s: %s" % (vtab,nam, sym.getinfo()))
#             
# 
#     def _show_groupstats(self,nam,st):
#         nam  = nam + ' '*(16-len(nam))
#         self.bprint('%s: %i variables, %i functions, %i subgroups' % (nam,st[0],st[1],st[2]))
# 
#     def _groupcount(g):
#         nvar, nfunc, ngroup = 0,0,0
#         for sym in g.values():
#             if  isGroup(sym):
#                 ngroup = ngroup + 1
#             elif sym.type in (symTypes.defpro,symTypes.pyfunc):
#                 nfunc  = nfunc  + 1
#             else: nvar = nvar + 1
#         return (nvar,nfunc, ngroup)
# 
