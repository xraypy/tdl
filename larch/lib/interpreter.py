'''Main Larch interpreter
'''
from __future__ import division, print_function
import os
import sys
import ast
try:
    import numpy
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False

import inputText
import builtins
from symboltable import SymbolTable, Group, isgroup
from util import Closure, LarchExceptionHolder, Procedure, DefinedVariable

__version__ = '0.9.2'
OPERATORS = {ast.Is:     lambda a, b: a is b,
             ast.IsNot:  lambda a, b: a is not b,
             ast.In:     lambda a, b: a in b,
             ast.NotIn:  lambda a, b: a not in b,
             ast.Add:    lambda a, b: a + b,
             ast.BitAnd: lambda a, b: a & b,
             ast.BitOr:  lambda a, b: a | b,
             ast.BitXor: lambda a, b: a ^ b,
             ast.Div:    lambda a, b: a / b,
             ast.FloorDiv: lambda a, b: a // b,
             ast.LShift: lambda a, b: a << b,
             ast.RShift: lambda a, b: a >> b,
             ast.Mult:   lambda a, b: a * b,
             ast.Pow:    lambda a, b: a ** b,
             ast.Sub:    lambda a, b: a - b,
             ast.Mod:    lambda a, b: a % b,
             ast.And:    lambda a, b: a and b,
             ast.Or:     lambda a, b: a or b,
             ast.Eq:     lambda a, b: a == b,
             ast.Gt:     lambda a, b: a > b,
             ast.GtE:    lambda a, b: a >= b,
             ast.Lt:     lambda a, b: a < b,
             ast.LtE:    lambda a, b: a <= b,
             ast.NotEq:  lambda a, b: a != b,
             ast.Invert: lambda a: ~a,
             ast.Not:    lambda a: not a,
             ast.UAdd:   lambda a: +a,
             ast.USub:   lambda a: -a}

class Interpreter:
    """larch program compiler and interpreter.
  This module compiles expressions and statements to AST representation,
  using python's ast module, and then executes the AST representation
  using a custom SymbolTable for named object (variable, functions).
  This then gives a restricted version of Python, with slightly modified
  namespace rules.  The program syntax here is expected to be valid Python,
  but that may have been translated as with the inputText module.

  The following Python syntax is not supported:
      Exec, Lambda, Class, Global, Generators, Yield, Decorators
        
  In addition, Function is greatly altered so as to allow a Larch procedure.
  """
    def __init__(self, symtable=None, writer=None):
        self.__writer = writer or sys.stdout.write
       
        if symtable is None:
            symtable = SymbolTable(larch=self)
        self.symtable   = symtable
        self._interrupt = None
        self.error      = [] 
        self.expr       = None
        self.retval     = None
        self.fname     = '<StdInput>'
        self.lineno    = -5
        builtingroup = getattr(symtable,'_builtin')
        mathgroup    = getattr(symtable,'_math')

        for sym in builtins.from_builtin:
            setattr(builtingroup, sym, __builtins__[sym])

        if HAS_NUMPY:
            for sym in builtins.from_numpy:
                setattr(mathgroup, sym, getattr(numpy, sym))
            
                for fname, sym in list(builtins.numpy_renames.items()):
                    setattr(mathgroup, fname, getattr(numpy, sym))

        for fname, fcn in list(builtins.local_funcs.items()):
            setattr(builtingroup, fname,
                    Closure(func=fcn, larch=self))
        setattr(builtingroup, 'definevar',
                Closure(func=self.__definevar))
        
    def __definevar(self, name, expr):
        """define a defined variable (re-evaluate on access)"""
        self.symtable.setSymbol(name, DefinedVariable(expr=expr,
                                                      larch=self))

    def unimplemented(self, node):
        "unimplemented nodes"
        cname = node.__class__.__name__
        self.on_except(node, "'%s' not supported" % (cname),
                       py_exc=sys.exc_info())

    def on_except(self, node, msg='', expr=None,
                  fname=None, lineno=-1, py_exc=None):
        "add an exception"
        if self.error is None:
            self.error = []
        if expr  is None:
            expr  = self.expr
        if fname is None:
            fname = self.fname        
        if lineno is None:
            lineno = self.lineno

        if len(self.error) > 0 and not isinstance(node, ast.Module):
            msg = 'Extra Error (%s)' % msg

        if py_exc is None:
            etype, evalue = None, None
        else:
            etype, evalue, tback = py_exc

        err = LarchExceptionHolder(node, msg=msg, expr= expr,
                                   fname= fname, lineno=lineno,
                                   py_exc=(etype, evalue) )
        self._interrupt = ast.Break()
        self.error.append(err)
    
    # main entry point for Ast node evaluation
    #  compile:  string statement -> ast
    #  interp :  ast -> result
    #  eval   :  string statement -> result = interp(compile(statement))
    def compile(self, text, fname=None, lineno=-4):
        """compile statement/expression to Ast representation    """
        self.expr  = text
        try:
            return  ast.parse(text)
        except:
            self.on_except(None, msg='Syntax Error',
                           expr=text, fname=fname, lineno=lineno,
                           py_exc=sys.exc_info())
            
    def interp(self, node, expr=None, fname=None, lineno=None):
        """executes compiled Ast representation for an expression"""
        # Note: keep the 'node is None' test: internal code here may run
        #    interp(None) and expect a None in return.
        if node is None:
            return None
        if isinstance(node, str):
            node = self.compile(node)

        method = "do_%s" % node.__class__.__name__.lower()
        if lineno is not None:
            self.lineno = lineno
        if fname  is not None:
            self.fname  = fname
        if expr   is not None:
            self.expr   = expr
        
        ret = None
        try:
            fcn = getattr(self, method)
        except:
            self.on_except(node, msg='Lookup Error', expr=expr,
                           fname=fname, lineno=lineno,
                              py_exc=sys.exc_info())
            return ret
        try:
            ret = fcn(node)
        except:
            self.on_except(node, msg='Runtime Error',
                              expr=expr, fname=fname, lineno=lineno,
                              py_exc=sys.exc_info())              
        if isinstance(ret, enumerate):
            ret = list(ret)
        return ret
            
    def eval(self, expr, fname=None, lineno=0):
        """evaluates a single statement"""
        self.fname = fname        
        self.lineno = lineno
        self.error = []
        node = self.compile(expr, fname=fname, lineno=lineno)
        if not self.error:
            out = self.interp(node, expr=expr,
                              fname=fname, lineno=lineno)
            if len(self.error) > 0:
                self.on_except(node, msg='Eval Error', expr=expr,
                               fname=fname, lineno=lineno,
                               py_exc=sys.exc_info())
            return out
        
    def dump(self, node, **kw):
        "simple ast dumper"
        return ast.dump(node, **kw)

    # handlers for ast components
    def do_expr(self, node):
        "expression"
        return self.interp(node.value)  # ('value',)

    def do_index(self, node):
        "index"
        return self.interp(node.value)  # ('value',)

    def do_return(self, node): # ('value',)
        "return statement"
        self.retval = self.interp(node.value)
        return
    
    def do_repr(self, node):
        "repr "
        return repr(self.interp(node.value))  # ('value',)

    def do_module(self, node):    # ():('body',) 
        "module def"
        out = None
        for tnode in node.body:
            out = self.interp(tnode)
        return out

    def do_expression(self, node):
        "basic expression"
        return self.do_module(node) # ():('body',) 

    def do_pass(self, node):
        "pass statement"
        return None  # () 

    def do_ellipsis(self, node):
        "ellipses"
        return Ellipsis

    # for break and continue: set the instance variable _interrupt
    def do_interrupt(self, node):    # ()
        "interrupt handler"
        self._interrupt = node
        return node

    def do_break(self, node):
        "break"
        return self.do_interrupt(node)

    def do_continue(self, node):
        "continue"
        return self.do_interrupt(node)    

    def do_assert(self, node):    # ('test', 'msg')
        "assert statement"
        if not self.interp(node.test):
            raise AssertionError(self.interp(node.msg()))
        return True

    def do_list(self, node):    # ('elt', 'ctx')
        "list"
        return [self.interp(e) for e in node.elts]

    def do_tuple(self, node):    # ('elts', 'ctx')
        "tuple"
        return tuple(self.do_list(node))
    
    def do_dict(self, node):    # ('keys', 'values')
        "dictionary"
        nodevals = list(zip(node.keys, node.values))
        interp = self.interp
        return dict([(interp(k), interp(v)) for k, v in nodevals])

    def do_num(self, node):
        'return number'
        return node.n  # ('n',) 

    def do_str(self, node):
        'return string'
        return node.s  # ('s',)

    def do_name(self, node):    # ('id', 'ctx')
        """ Name node """
        ctx = node.ctx.__class__
        if ctx == ast.Del:
            val = self.symtable.delSymbol(node.id)
        elif ctx == ast.Param:  # for Function Def
            val = str(node.id)
        else:
            val = self.symtable.getSymbol(node.id)
            if isinstance(val, DefinedVariable):
                val = val.evaluate()
        return val

    def node_assign(self, nod, val):
        """here we assign a value (not the node.value object) to a node
        this is used by do_assign, but also by for, list comprehension, etc.
        """
        if len(self.error) > 0:
            return
        if nod.__class__ == ast.Name:
            sym = self.symtable.setSymbol(nod.id, value=val)
        elif nod.__class__ == ast.Attribute:
            if nod.ctx.__class__  == ast.Load:
                errmsg = "cannot assign to attribute %s" % nod.attr
                self.on_except(nod, errmsg)
            setattr(self.interp(nod.value), nod.attr, val)
        elif nod.__class__ == ast.Subscript:
            sym    = self.interp(nod.value)
            xslice = self.interp(nod.slice)
            if isinstance(nod.slice, ast.Index):
                sym.__setitem__(xslice, val)
            elif isinstance(nod.slice, ast.Slice):
                sym.__setslice__(xslice.start, xslice.stop, val)
            elif isinstance(nod.slice, ast.ExtSlice):
                sym[(xslice)] = val
        elif nod.__class__ in (ast.Tuple, ast.List):
            if len(val) == len(nod.elts):
                for telem, tval in zip(nod.elts, val):
                    self.node_assign(telem, tval)
            else:
                raise ValueError('too many values to unpack')

    def do_attribute(self, node):    # ('value', 'attr', 'ctx')
        "extract attribute"
        ctx = node.ctx.__class__
        # print("do_attribute",node.value,node.attr,ctx)
        if ctx == ast.Load:
            sym = self.interp(node.value)
            if hasattr(sym, node.attr):
                val = getattr(sym, node.attr)
                if isinstance(val, DefinedVariable):
                    val = val.evaluate()
                return val
            else:
                obj = self.interp(node.value)
                fmt = "%s does not have member '%s'"                
                if not isgroup(obj):
                    obj = obj.__class__
                    fmt = "%s does not have attribute '%s'"
                msg = fmt % (obj, node.attr)
                self.on_except(node, msg=msg, py_exc=sys.exc_info())
        elif ctx == ast.Del:
            return delattr(sym, node.attr)
        elif ctx == ast.Store:
            msg = "attribute for storage: shouldn't be here!"
            self.on_except(node, msg=msg, py_exc=sys.exc_info())        

    def do_assign(self, node):    # ('targets', 'value')
        "simple assignment"
        val = self.interp(node.value)
        if len(self.error) > 0:
            return        
        for tnode in node.targets:
            self.node_assign(tnode, val)
        return # return val

    def do_augassign(self, node):    # ('target', 'op', 'value')
        "augmented assign"
        # print( "AugASSIGN ", node.target, node.value)
        return self.do_assign(ast.Assign(targets=[node.target],
                                         value=ast.BinOp(left = node.target,
                                                         op   = node.op,
                                                         right= node.value)))
       
    def do_slice(self, node):    # ():('lower', 'upper', 'step')
        "simple slice"
        return slice(self.interp(node.lower), self.interp(node.upper),
                     self.interp(node.step))

    def do_extslice(self, node):    # ():('dims',)
        "extended slice"
        return tuple([self.interp(tnode) for tnode in node.dims])
    
    def do_subscript(self, node):    # ('value', 'slice', 'ctx') 
        "subscript handling -- one of the tricky parts"
        # print("do_subscript: ", ast.dump(node))
        val    = self.interp(node.value)
        nslice = self.interp(node.slice)
        ctx = node.ctx.__class__
        if ctx in ( ast.Load, ast.Store):
            if isinstance(node.slice, (ast.Index, ast.Slice, ast.Ellipsis)):
                return val.__getitem__(nslice)
            elif isinstance(node.slice, ast.ExtSlice):
                return val[(nslice)]
        else:
            msg = "subscript with unknown context"
            self.on_except(node, msg=msg, py_exc=sys.exc_info())

    def do_delete(self, node):    # ('targets',)
        "delete statement"
        error = False
        for tnode in node.targets:
            ctx = tnode.ctx.__class__
            if ctx != ast.Del:
                error = True
                break

            children = []
            while tnode.__class__ == ast.Attribute:
                children.append(tnode.attr)
                tnode = tnode.value

            if tnode.__class__ == ast.Name:
                parent,child = self.symtable.parentOf(tnode.id)
                children.append(child)
                lastchild = children.pop(0)
                for child in children:
                    if hasattr(parent, child):
                        parent = getattr(parent, child)
                    else:
                        error = True
                        break
                delattr(parent, lastchild)
            else:
                error = True
                break
        if error:
            msg = "could not delete symbol"
            self.on_except(node, msg=msg, py_exc=sys.exc_info())
                
            
    def do_unaryop(self, node):    # ('op', 'operand')
        "unary operator"
        return OPERATORS[node.op.__class__](self.interp(node.operand))
    
    def do_binop(self, node):    # ('left', 'op', 'right')
        "binary operator"
        return OPERATORS[node.op.__class__](self.interp(node.left),
                                            self.interp(node.right))

    def do_boolop(self, node):    # ('op', 'values')
        "boolean operator"
        val = self.interp(node.values.pop(0))
        is_and = ast.Or != node.op.__class__
        if (is_and and val) or (not is_and and not val):
            for n in node.values:
                val =  OPERATORS[node.op.__class__](val, self.interp(n))
                if (is_and and not val) or (not is_and and val):
                    break
        return val
    
    def do_compare(self, node):    # ('left', 'ops', 'comparators')
        "comparison operators"
        lval = self.interp(node.left)
        out  = True
        for oper, rnode in zip(node.ops, node.comparators):
            comp = OPERATORS[oper.__class__]
            rval = self.interp(rnode)
            out  = out and  comp(lval, rval)
            lval = rval
            if not out:
                break
        return out

    def do_print(self, node):    # ('dest', 'values', 'nl')
        """ note: implements Python2 style print statement, not
        print() function.  Probably, the 'larch2py' translation
        should look for and translate print -> print_() to become
        a customized function call.
        """
        dest = self.interp(node.dest) or sys.stdout
        end = ''
        if node.nl:
            end = '\n'
        out = [self.interp(tnode) for tnode in node.values]
        if out and len(self.error)==0:
            print(*out, file=dest, end=end)
        
    def do_if(self, node):    # ('test', 'body', 'orelse')
        "regular if-then-else statement"
        block = node.orelse
        if self.interp(node.test):
            block = node.body
        for tnode in block:
            self.interp(tnode)

    def do_ifexp(self, node):    # ('test', 'body', 'orelse')
        "if expressions"
        expr = node.orelse
        if self.interp(node.test):
            expr = node.body
        return self.interp(expr)

    def do_while(self, node):    # ('test', 'body', 'orelse')
        "while blocks"
        while self.interp(node.test):
            self._interrupt = None
            for tnode in node.body:
                self.interp(tnode)
                if self._interrupt is not None:
                    break
            if isinstance(self._interrupt, ast.Break):
                break
        else:
            for tnode in node.orelse:
                self.interp(tnode)
        self._interrupt = None

    def do_for(self, node):    # ('target', 'iter', 'body', 'orelse')
        "for blocks"
        for val in self.interp(node.iter):
            self.node_assign(node.target, val)
            if len(self.error) > 0:
                return            
            self._interrupt = None
            for tnode in node.body:
                self.interp(tnode)
                if len(self.error) > 0:
                    return                
                if self._interrupt is not None:
                    break
            if isinstance(self._interrupt, ast.Break):
                break
        else:
            for tnode in node.orelse:
                self.interp(tnode)
        self._interrupt = None

    def do_listcomp(self, node):    # ('elt', 'generators') 
        "list comprehension"
        out = []
        for tnode in node.generators:
            if tnode.__class__ == ast.comprehension:
                for val in self.interp(tnode.iter):
                    self.node_assign(tnode.target, val)
                    if len(self.error) > 0:
                        return                    
                    add = True
                    for cond in tnode.ifs:
                        add = add and self.interp(cond)
                    if add:
                        out.append(self.interp(node.elt))
        return out

    def do_call(self, node):
        "function/procedure execution"
        # ('func', 'args', 'keywords', 'starargs', 'kwargs')
        func = self.interp(node.func)
        #  note: callable(func) in 2.x becomes
        #  hasattr(func, '__call__') in 3.x, BUT::
        #  classes in 2.x do NOT HAVE attribute
        # '__call__'  and ARE callable!!
        if not callable(func): # hasattr(func, '__call__'):
            msg = "'%s' is not callable" % (func)
            self.on_except(node, msg=msg, py_exc=sys.exc_info())

        args = [self.interp(targ) for targ in node.args]
        if node.starargs is not None:
            args = args + self.interp(node.starargs)
        
        keywords = {}
        for key in node.keywords:
            if not isinstance(key, ast.keyword):
                msg = "keyword error in function call '%s'" % (func)
                self.on_except(node, msg=msg, py_exc=sys.exc_info())
            
            keywords[key.arg] = self.interp(key.value)
        if node.kwargs is not None:
            keywords.update(self.interp(node.kwargs))
        return func(*args, **keywords)
    
    def do_functiondef(self, node):
        "define procedures"
        # ('name', 'args', 'body', 'decorator_list') 
        if node.decorator_list != []:
            print("Warning: decorated procedures not supported!")

        kwargs = []
        while node.args.defaults:
            defval = self.interp(node.args.defaults.pop())
            key    = self.interp(node.args.args.pop())
            kwargs.append((key, defval))
        kwargs.reverse()
        args = [tnode.id for tnode in node.args.args]
        doc = None
        if isinstance(node.body[0], ast.Expr):
            docnode = node.body.pop(0)
            doc = self.interp(docnode.value)
        # 
        proc = Procedure(node.name, larch= self, doc= doc,
                         body   = node.body,
                         fname  = self.fname,   lineno = self.lineno,
                         args   = args,   kwargs = kwargs,
                         vararg = node.args.vararg,
                         varkws = node.args.kwarg)
        self.symtable.setSymbol(node.name, value=proc)

    # imports
    def do_import(self, node):    # ('names',)
        "simple import"
        for tnode in node.names:
            self.import_module(tnode.name, asname=tnode.asname)
        
    def do_importfrom(self, node):    # ('module', 'names', 'level')
        "import/from"
        fromlist, asname = [], []
        for tnode in node.names:
            fromlist.append(tnode.name)
            asname.append(tnode.asname)
        self.import_module(node.module,
                           asname=asname, fromlist=fromlist)


    def import_module(self, name, asname=None,
                      fromlist=None, do_reload=False):
        """
        import a module (larch or python), installing it into the symbol table.
        required arg:
            name       name of module to import
                          'foo' in 'import foo'
        options:
            fromlist   list of symbols to import with 'from-import'
                          ['x','y'] in 'from foo import x, y'
            asname     alias for imported name(s)
                          'bar' in 'import foo as bar'
                       or
                          ['s','t'] in 'from foo import x as s, y as t'

        this method covers a lot of cases (larch or python, import
        or from-import, use of asname) and so is fairly long.
        """
        symtable = self.symtable
        st_sys     = symtable._sys
        
        for idir in st_sys.path:
            if idir not in sys.path and os.path.exists(idir):
                sys.path.append(idir)

        # step 1  import the module to a global location
        #   either sys.modules for python modules
        #   or  st_sys.modules for larch modules
        # reload takes effect here in the normal python way:
        #   if
        do_load = (name not in st_sys.modules and name not in sys.modules)
        do_load = do_load or do_reload
        # print("  larch do_load ", do_load, len(self.error))
        if do_load:
            # first look for "name.lar"
            islarch = False
            larchname = "%s.lar" % name
            for dirname in st_sys.path:
                if not os.path.exists(dirname):
                    continue
                if larchname in os.listdir(dirname):
                    islarch = True
                    modname = os.path.abspath(os.path.join(dirname, larchname))
                    # print("Found larch module:" , modname)
                    # save current module group
                    #  create new group, set as moduleGroup and localGroup
                    symtable.save_frame()
                    st_sys.modules[name] = thismod = Group(name=name)
                    symtable.set_frame((thismod, thismod))

                    text = open(modname).read()
                    inptext = inputText.InputText()
                    inptext.put(text, filename=modname)

                    while inptext:
                        block, fname, lineno = inptext.get()
                        self.eval(block, fname=fname, lineno=lineno)
                        if self.error:
                            break
                    symtable.restore_frame()
            if self.error:
                st_sys.modules.pop(name)
                thismod = None
                return

            # or, if not a larch module, load as a regular python module
            if not islarch:
                try:
                    __import__(name)
                    thismod = sys.modules[name]
                except:
                    # print(" import error ", name, sys.exc_info())
                    self.on_except(None, msg='Import Error',
                                      py_exc=sys.exc_info())
                    return
        else: # previously loaded module, just do lookup
            if name in st_sys.modules:
                thismod = st_sys.modules[name]
            elif name in sys.modules:
                thismod = sys.modules[name]               
               
        # now we install thismodule into the current moduleGroup
        # import full module
        if fromlist is None:
            if asname is None:
                asname = name
            parts = asname.split('.')
            asname = parts.pop()
            top = st_sys.moduleGroup
            while len(parts) > 0:
                subname = parts.pop(0)
                subgrp  = Group()
                setattr(top, subname, subgrp)
                top = subgrp
            setattr(top, asname, thismod)
        # import-from construct
        else:
            if asname is None:
                asname = [None]*len(fromlist)
            for sym, alias in zip(fromlist, asname):
                if alias is None:
                    alias = sym
                setattr(st_sys.moduleGroup, alias, getattr(thismod, sym))
    # end of import_module

    # not yet implemented:
    def do_excepthandler(self, node): # ('type', 'name', 'body')
        "exception handler..."
        print("except handler %s / %s " % (node.type, ast.dump(node.name)))
        return (self.interp(node.type), node.name, node.body)
    
    def do_tryexcept(self, node):    # ('body', 'handlers', 'orelse')
        "try/except blocks"
        for tnode in node.body:
            self.interp(tnode)
            if self.error:
                e_type, e_value = self.error[-1].py_exc
                this_exc = e_type()
                # print("Look for except: ", this_exc)
                # print("out of handlers: ", node.handlers)
                for hnd in node.handlers:
                    htype = None
                    if hnd.type is not None:
                        htype = __builtins__.get(hnd.type.id, None)
                    if htype is None or isinstance(this_exc, htype):
                        self.error = []
                        if hnd.name is not None:
                            self.node_assign(hnd.name, e_value)
                        for tline in hnd.body:
                            self.interp(tline)
                        break

    def do_raise(self, node):    # ('type', 'inst', 'tback')
        "raise statement"
        msg = "%s: %s" % (self.interp(node.type).__name__,
                          self.interp(node.inst))
        self.on_except(node.type, msg=msg,
                       py_exc=sys.exc_info())
                    
    def do_tryfinally(self, node):
        "try finally"
        return self.unimplemented(node)
    
    def do_exec(self, node):
        "exec statement: not allowed!"
        return self.unimplemented(node)
    
    def do_lambda(self, node):
        "lambda: non-ultimate"
        return self.unimplemented(node)
    
    def do_class(self, node):
        "class: we have none!"
        return self.unimplemented(node)
    
    def do_global(self, node):
        "global: handled with different namespace rules"
        return self.unimplemented(node)
    
    def do_generators(self, node):
        "generators: not implemented"
        return self.unimplemented(node)

    def do_yield(self, node):
        "yield: not implemented"
        return self.unimplemented(node)

    def do_decorators(self, node):
        "decorators: not implemented"
        return self.unimplemented(node)
    
    def do_generatorexp(self, node):
        "generator expressions: not implemented"
        return self.unimplemented(node)
    
