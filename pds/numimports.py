#################################################################
"""
Tom Trainor (fftpt@uaf.edu), 11-08
Numerical/plotting features for shell program

Modifications:
--------------
- Initial ideas from old tdl NumImports module by Matt Newville

"""
#################################################################
"""
Notes:
The _NumShell class should be used as a subclass within shell.Shell.
This adds numerical/plotting modules to the shell by default,
as well as the calc function.

If numpy is not installed this should still be ok to use as a subclass
of shell.Shell, it just wont add any of the numerical stuff... 

The numerical/plotting modules are handled by the function call:
(num,scipy,pylab) = _modules(backend="Wx")

"""
#################################################################

import sys
import types
from   shellutil import mod_import, show_more, PrintExceptErr

#################################################################

#################################################################
# define some utility/convience functions
#################################################################

#################################################################
def _parse_version(s):
    """ parse a version string to an integer '1.0.1' -> 101"""
    factor = 100.
    version = 0.0
    vlst = []
    #for v in [int(i) for i in s.split('.')]:
    for v in [i for i in s.split('.')]:
        try:
            version = version + int(v) * factor
        except:
            pass
        factor = factor * 0.1000
    return version

#################################################################
def _modules(backend="TkAgg"):
    """
    import numpy, scipy and pylab - return as a tuple.
    (num,scipy,pylab) = Modules(backend="TkAgg")
    """
    ## numpy
    num_version = 0
    try:
        num         = mod_import('numpy')
        version     = _parse_version(num.__version__)
        num_version = 'numpy %s' % num.__version__
        # add some stuff
        num.ArrayType = num.ndarray
    except:
        #raise ImportError, 'Error importing numpy'
        s = "Warning no numpy"
        PrintExceptErr(s)
        num = None
    if num_version < 100:
        raise ImportError, 'Need numpy version 1.0 or higher: %s' % num_version    

    ## scipy 
    try:
        #import scipy
        scipy   = mod_import('scipy')
        version = _parse_version(scipy.__version__)
        scipy_version = 'scipy %s ' % version
        if version < 48:
            print "Warning Scipy version is old: %s" % scipy_version
    except:
        print "Warning no scipy"
        scipy = None
 
    ## pylab 
    try:
        pylab = _import_pylab(backend=backend)
    except:
        s = "Warning no pylab"
        PrintExceptErr(s)
        pylab = None

    return (num,scipy,pylab)

#################################################################
def _import_pylab(backend="TkAgg",verbose=True):
    #
    import matplotlib
    check = matplotlib.get_backend()
    #
    if backend == "TkAgg":
        # below is kluge for compiled version, ie the way we are 
        # running py2exe -> basically forces WX as backend.
        #if check == "WXAgg":
        #    pass
        #else:
        import Tkinter
        PLOT_ROOT = Tkinter.Tk()
        PLOT_ROOT.iconify()
    elif backend == "WXAgg":
        import wx
        PLOT_ROOT = None
    #
    version = _parse_version(matplotlib.__version__)
    txt = "    ** INIT MATPLOTLIB, backend = %s, version=%3.1f\n\n" % (backend,version)
    if verbose: sys.__stdout__.write(txt)
    #
    if version < 950:
        matplotlib.use(backend)
    else:
        matplotlib.use(backend,warn=False)
    import matplotlib.pylab as pylab
    
    # do we need this??
    # ie plotting works without below line
    # but this lets us grab tk if needed...
    pylab.plot_root = PLOT_ROOT 
    
    # import plotter
    import plotter
    pylab.plotter    = plotter.plotter
    pylab.newplotter = plotter.newplotter
    pylab.cursor     = plotter.cursor

    # Make sure we're interactive
    #pylab.show._needmain=False
    #matplotlib.interactive(True)
    pylab.ion()
    
    return pylab

#################################################################
#################################################################
class _NumShell:
    
    def num_setup(self,):
        """
        Add numerical methods to Shell.shell
        """
        ##############
        self.pds_commands.update({'calc':self.do_calc})

        #############        
        # add modules to symbol table
        (num,scipy,pylab) = _modules(backend=self.GUI)

        if num: 
            self.interp.symbol_table.data["num"] = num

        if scipy:
            self.interp.symbol_table.data["scipy"] = scipy
            #self.interp.symbol_table.data["optimize"] = __import__("scipy.optimize")

        if pylab:
            self.interp.symbol_table.data["pylab"] = pylab
            #self.close_pylab = pylab.close
            def _close_pylab():
                pylab.close('all')
            self.close_pylab = _close_pylab

        # add to startup commands
        startup = []

        # put some common numpy stuff into __builtins__
        if num:
            startup.append("__builtins__.update({'sin':num.sin})")
            startup.append("__builtins__.update({'asin':num.arcsin})")
            startup.append("__builtins__.update({'cos':num.cos})")
            startup.append("__builtins__.update({'acos':num.arccos})")
            startup.append("__builtins__.update({'tan':num.tan})")
            startup.append("__builtins__.update({'atan':num.arctan})")
            startup.append("__builtins__.update({'deg':num.degrees})")
            startup.append("__builtins__.update({'rad':num.radians})")
            
            startup.append("__builtins__.update({'ln':num.log})")
            startup.append("__builtins__.update({'log':num.log10})")
            startup.append("__builtins__.update({'exp':num.exp})")
            startup.append("__builtins__.update({'sqrt':num.sqrt})")
            startup.append("__builtins__.update({'pi':num.pi})")
            startup.append("__builtins__.update({'e':num.e})")

        # put some mathutil functions into __builtins__
        try:
            import mpcutils.mathutil as __mathutil__
            self.interp.symbol_table.data["__mathutil__"] = __mathutil__
            startup.append("__builtins__.update({'ave':__mathutil__.ave})")
            startup.append("__builtins__.update({'std':__mathutil__.std})")
            startup.append("__builtins__.update({'square':__mathutil__.square})")
            startup.append("__builtins__.update({'minimize':__mathutil__.minimize})")
            #
            startup.append("__builtins__.update({'cosd':__mathutil__.cosd})")
            startup.append("__builtins__.update({'acosd':__mathutil__.arccosd})")
            startup.append("__builtins__.update({'sind':__mathutil__.sind})")
            startup.append("__builtins__.update({'asind':__mathutil__.arcsind})")
            startup.append("__builtins__.update({'tand':__mathutil__.tand})")
            startup.append("__builtins__.update({'atand':__mathutil__.arctand})")
            #
            startup.append("__builtins__.update({'random':__mathutil__.random})")
            startup.append("__builtins__.update({'voigt':__mathutil__.voigt})")
            #
            startup.append("__builtins__.update({'Peak':__mathutil__.Peak})")
            startup.append("__builtins__.update({'LinReg':__mathutil__.LinReg})")
        except:
            print "Cant load mathutil functions"
        
        # put some plotter stuff into __builtins__
        if pylab:
            startup.append("__builtins__.update({'plot':pylab.plotter})")
            startup.append("__builtins__.update({'newplot':pylab.newplotter})")
            startup.append("addcmd 'plot', 'pylab.plotter'")
            startup.append("addcmd 'newplot', 'pylab.newplotter'")

        return startup

    ##############################################################

    #############################################################
    #  Calculator
    #############################################################
        
    ##############################################################
    def do_calc(self, line=''):
        """
        Calculator.  Help is available from the 'calc' command prompt
        """
        help_1 = "*** Calculator: 'q' = exit, 'h' for help, 'u' = use examples"
        help_2 = """*** Calculator Commands
        \r'c' = clear buffer
        \r'd' = debug (default off)
        \r'f' = convert results to float (default off)
        \r'h' = help
        \r'p' = pop last entry from buffer
        \r'q' = exit
        \r's' = show buffer
        \r'u' = useage examples
        \r'w' = swap last two entries in buffer
        """
        use = """
        \r*** Useage
        \r In the below examples the result of a computation is
        \r appended to the buffer.  Use:
        \r   's' to inspect the buffer,
        \r   'w' to swap the top 2 buffer values,
        \r   'p' to pop (remove) the top buffer value
        \r Note, be cautious of integer arithmetic.  Even if the
        \r float flag is turned on ('f'), only the result
        \r of an evaluated expression will be converted to float. 
        \r Therefore, its best to explicitly use floating point numbers...
        
        \r*** Examples
        \r# Place values into the buffer (note b is the top of the buffer)
        \rcalc>3                # --> a
        \rcalc>6                # --> b  
        
        \r# A sole operate works on last two buffer entries
        \rcalc>-                # --> op
        \r  = -3                #  3 - 6 = a op b

        \r# A trailing operator works the same
        \rcalc>3                # --> a
        \rcalc>6-               # --> b op
        \r  = -3                #  3 - (6) = a op (b)

        \r# When an operator trails an expression the buffer
        \r# the expression is evaluated before the opertor
        \r# and buffer value are combined, i.e.
        \r# the expression is placed in parenthesis
        \rcalc>3                # --> a
        \rcalc>6+4-             # --> expr op
        \r  = -7                #  3 - (6 + 4) = a op (expr) 

        \r# The top value in the buffer can be used in an
        \r# expression with pop:
        \rcalc>10               # --> a
        \rcalc>num.log10(pop)   # --> expr(a)
        \r  = 1.0               #  num.log(3) = expr(a)

        \r# The pop occurs first, so when combined with the
        \r# operator syntax, the buffer value in the constructed
        \r# expression is the second in the list (note only one
        \r# pop is allowed in an expression)
        \rcalc>1                # --> a
        \rcalc>10               # --> b
        \rcalc>num.log10(pop)-  # --> expr(b) op
        \r  = 0.0               # 1 - (num.log10(10)) = a op (expr(b))
        \r
        """
        ################
        ex = self.interp.execute
        ev = self.interp.evaluate
        
        ################
        if len(line) > 0:
            print ev(line)
            return

        ################
        print help_1
        self.interp.symbol_table.del_symbol('__calc__')
        ex('__calc__ = group()')
        self.interp.set_data('__calc__.buff',[])
        self.interp.set_data('__calc__.val',0.0)

        ################
        def buff_len():
            tmp = self.interp.get_data('__calc__.buff')
            if tmp == None: return 0
            return len(tmp)
        
        ################
        def get_line(line):
            # check for pop's
            idx = line.find('pop')
            if idx > -1:
                if buff_len() > 0:
                    ex("__calc__.p = __calc__.buff.pop()",print_err=False)
                    line = line[:idx] + "__calc__.p" + line[idx+3:]
                else:
                    print "Buffer empty"
                    return None

            # look at operators...
            if line in ('+','-','/','*'):
                """
                here do: former op later --> a op b
                """
                if buff_len() < 2:
                    print "Buffer has fewer than 2 values"
                    line = None
                else:
                    ex("__calc__.b = __calc__.buff.pop()",print_err=False)
                    ex("__calc__.a = __calc__.buff.pop()",print_err=False)
                    line = "__calc__.a  %s  __calc__.b" % (line)
            elif line[-1] in ('+','-','/','*'):
                """
                here do: buff op expr 
                """
                op = line[-1]
                line = line[:-1]
                if buff_len() < 1:
                    print "Buffer empty"
                    line = None
                else:
                    ex("__calc__.a = __calc__.buff.pop()",print_err=False)
                    line = "__calc__.a %s (%s)" % (op,line)
            return line
        
        ################
        debug = False
        do_float = False
        calc_types = [types.BooleanType, types.ComplexType,
                      types.FloatType, types.IntType,
                      types.LongType]
        try:
            import numpy
            calc_types.append([numpy.ndarray])
            calc_types.append(numpy.typeDict.values())
        except:
            numpy = None
            
        ################
        while 1:
            line = raw_input('calc>')
            line = line.strip()
            if len(line) > 0:
                if line in ('q','quit'):
                    return
                elif line in ('h','help'):
                    print help_2
                elif line in ('u','use'):
                    #print use
                    show_more(use)
                elif line in ('d','debug'):
                    debug = not debug
                    print "Debug = ", debug
                elif line in ('f','float'):
                    do_float = not do_float
                    print "Do float = ", do_float
                elif line in ('c','clear'):
                    ex('__calc__.buff = []')
                elif line in ('s','show'):
                    buff = self.interp.get_data('__calc__.buff')
                    for v in buff: print v
                elif line in ('w','swap'):
                    if buff_len() > 1:
                        ex('__calc__.a = __calc__.buff.pop()')
                        ex('__calc__.b = __calc__.buff.pop()')
                        ex('__calc__.buff.append(__calc__.a)')
                        ex('__calc__.buff.append(__calc__.b)')
                elif line in ('p','pop'):
                    if buff_len() > 0:
                        ex('print __calc__.buff.pop()')
                else:
                    # finally execute line
                    line = get_line(line)
                    if debug: print line
                    if line:
                        if debug: print_err = True
                        else: print_err = False
                        val = ev(line,print_err=print_err)
                        if type(val) in calc_types:
                            if do_float:
                                try:
                                    val = float(val)
                                except:
                                    pass
                            print "  = " ,  val
                            self.interp.set_data('__calc__.val',val)
                            ex('__calc__.buff.append(__calc__.val)')
                        else:
                            print "Error in type returned to calc: %s" % type(val)
                    else:
                        val = None
                        self.interp.set_data('__calc__.val',0.0)

#################################################################
