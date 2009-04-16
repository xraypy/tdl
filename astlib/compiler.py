from __future__ import division, print_function
import os
import sys
import ast
import inputText
import symbolTable
import builtins
import numpy
import __builtin__
from util import closure

class TDLError(Exception):
    def __init__(self,error,descr = None,node = None):
        self.error = error
        self.descr = descr
    def __repr__(self):
        return "%s: %s" % (self.error, self.descr)
    __str__ = __repr__

class TDLStopError(Exception):
    def __init__(self,error,descr = None,node = None):
        self.error = error
        self.descr = descr
    def __repr__(self):
        return "%s: %s" % (self.error, self.descr)
    __str__ = __repr__


__version__ = '0.9.0'

major,minor = sys.version_info[0], sys.version_info[1]
if major == 2 and minor < 6:
    raise EnvironmentError('requires python 2.6 or higher')

_operators = {ast.Is:     lambda a,b: a is b,
              ast.IsNot:  lambda a,b: a is not b,
              ast.In:     lambda a,b: a in b,
              ast.NotIn:  lambda a,b: a not in b,
              ast.Add:    lambda a,b: a + b,
              ast.BitAnd: lambda a,b: a & b,
              ast.BitOr:  lambda a,b: a | b,
              ast.BitXor: lambda a,b: a ^ b,
              ast.Div:    lambda a,b: a / b,
              ast.FloorDiv: lambda a,b: a // b,
              ast.LShift: lambda a,b: a << b,
              ast.RShift: lambda a,b: a >> b,
              ast.Mult:   lambda a,b: a * b,
              ast.Pow:    lambda a,b: a ** b,
              ast.Sub:    lambda a,b: a - b,
              ast.Mod:    lambda a,b: a % b,
              ast.And:    lambda a,b: a and b,
              ast.Or:     lambda a,b: a or b,
              ast.Eq:     lambda a,b: a==b,
              ast.Gt:     lambda a,b: a>b,
              ast.GtE:    lambda a,b: a>=b,
              ast.Lt:     lambda a,b: a<b,
              ast.LtE:    lambda a,b: a<=b,
              ast.NotEq:  lambda a,b: a!=b,
              ast.Invert: lambda a: ~a,
              ast.Not:    lambda a: not a,
              ast.UAdd:   lambda a: +a,
              ast.USub:   lambda a: -a }

def _Class(node):     return node.__class__
def _ClassName(node): return node.__class__.__name__
def _Context(node):   return node.ctx.__class__
def _Op(node):        return _operators[node.__class__]

##
class DefinedVariable(object):
    """defined variable: re-evaluate on access

    Note that the localGroup/moduleGroup are cached
    at compile time, and restored for evaluation.
    """
    def __init__(self, expr=None, compiler=None):
        self.expr = expr
        self.compiler = compiler
        self.ast = None
        self._groups = None,None
        self.compile()

    def compile(self):
        if self.compiler is not None and self.expr is not None:
            self.ast = self.compiler.compile(self.expr)
            _sys = self.compiler.symtable._sys
            self._groups = (_sys.localGroup,_sys.moduleGroup)

    def evaluate(self):
        if self.ast is None: self.compile()

        if self.ast is None:
            msg = "Cannot compile '%s'"  % (self.expr)
            raise Warning, msg            
            
        if hasattr(self.compiler,'interp'):
            _sys = self.compiler.symtable._sys
            # save current localGroup/moduleGroup 
            save_groups  = _sys.localGroup,_sys.moduleGroup
            
            _sys.localGroup,_sys.moduleGroup = self._groups
            rval = self.compiler.interp(self.ast)

            _sys.localGroup,_sys.moduleGroup = save_groups
            return rval
        else:
            msg = "Cannot evaluate '%s'"  % (self.expr)
            raise ValueError, msg

def definevar(name,expr,compiler=None):
    if compiler is not None:
        defvar = DefinedVariable(expr=expr,compiler=compiler)
        compiler.symtable.setSymbol(name,defvar)
     
class Procedure(object):
    """tdl procedure:  function """
    def __init__(self, name, compiler=None,
                 body=None, args=None, kwargs=None,
                 vararg=None, varkws=None):
        self.name     = name
        self.compiler = compiler
        self.modgroup = compiler.symtable._sys.moduleGroup
        self.body     = body
        self.argnames = args
        self.kwargs   = kwargs
        self.vararg   = vararg
        self.varkws   = varkws
        
    def __call__(self,*args,**kwargs):
        stable  = self.compiler.symtable
        sys     = stable._sys      
        lgroup  = stable.createGroup()

        args = list(args)
        for argname in self.argnames:
            setattr(lgroup, argname,args.pop(0))
        if self.vararg is not None and len(args)>0:
            setattr(lgroup, self.vararg, args)
        
        for k,v in self.kwargs.items():
            if kwargs.has_key(k):  v = kwargs.pop(k)
            setattr(lgroup, k, v)

        if self.varkws is not None and len(kwargs)>0:
            setattr(lgroup, self.varkws, kwargs)
        
        grps_save = sys.localGroup,sys.moduleGroup
        stable._set_local_mod((lgroup, self.modgroup))
        
        retval = None
        self.compiler.retval = None
        for node in self.body:
            self.compiler.interp(node)
            if self.compiler.retval is not None:
                retval = self.compiler.retval
                break

        sys.localGroup,sys.moduleGroup = grps_save
        stable._set_local_mod(grps_save)
        self.compiler.retval = None
        del lgroup
        return retval
    
class ExceptionStorage:
    def __init__(self):
        self.raised = False
        self.node   = None
        self.msg    = []
        self.lineno = None
        self.col_offset = 0

        
class Compiler:
    """ program compiler and interpreter.
  This module compiles expressions and statements to AST representation,
  using python's ast module, and then executes the AST representation
  using a custom SymboplTable for named object (variable, functions).
  This then gives a restricted version of Python,
    
  The following Python syntax is not supported:
      Generators, Yield, Decorators, Class, Exec, Lambda, Global
        
  In addition, Function is greatly altered so as to allow 
    
    """
    def __init__(self,symtable=None,input=None, writer=None,
                 load_builtins=True):
        self.__writer = writer or sys.stdout.write
        if symtable is None:
            symtable = symbolTable.symbolTable(writer=self.__writer)
        self.symtable    = symtable
        self.setSymbol   = symtable.setSymbol
        self.getSymbol   = symtable.getSymbol
        self.delSymbol   = symtable.delSymbol        
        self._interrupt  = None
        self._exception_raised = False        
        self.error_msg   = []
        self.retval      = None

        imports = ((builtins._from_builtin, __builtin__,'_builtin'),
                   (builtins._from_numpy,   numpy, '_math'))
        
        for symlist, pymod, tdlmod in imports:
            group = getattr(symtable,tdlmod)
            for sym in symlist:
                setattr(group, sym, getattr(pymod, sym))

        group = getattr(symtable, '_builtin')
        for fname,fcn in builtins._local_funcs.items():
            setattr(group, fname, closure(func=fcn,compiler=self))
        setattr(group, 'definevar', closure(func=definevar,compiler=self))
        

    ##
    def NotImplemented(self,node):
        self.onError(TDLError, node, msg="'%s' not supported" % (_ClassName(node)))

    # main entry point for Ast node evaluation
    #  compile:  string statement -> ast
    #  interp :  ast -> result
    #  eval   :  string statement -> result = interp(compile(statement))
    def compile(self,text):
        """compile statement/expression to Ast representation    """
        try:
            return ast.parse(text)
        except:
            self.onError(SyntaxError, node)
            raise
        
    def interp(self, node):
        """executes compiled Ast representation for an expression"""

        # Note: keep the 'node is None' test: internal code here may run
        #    interp(None) and expect a None in return.
        if node is None: return None
        if isinstance(node,(str,unicode)):  node = self.compile(node)
        try:
            return getattr(self,"do%s" % _ClassName(node))(node)
        except:
            self.onError(TDLError, node)
            raise

            
    def eval(self,expr,fname=None,lineno=None):
        """evaluates a single statement"""
        self.fname = fname
        self.lineno = lineno
        self.error = ExceptionStorage()

        try:
            node = self.compile(expr)
        except:
            self.onError(TDLError, None)
            self.show_error(expr=expr)
        try:
            node = self.compile(expr)
            return self.interp(node)
        except:
            self.onError(TDLError, node)
            self.show_error(expr=expr)

    def show_error(self,expr=None):
        if expr is not None:
            elines = expr.split('\n')
            nlines = len(elines)
            if nlines > 1 and nlines > self.error.lineno:
                expr = elines[self.error.lineno-1]
        self.error.msg.append("%s" % expr)
        if self.error.col_offset >0:
            self.error.msg.append('%s^' % (' '*self.error.col_offset))
        while self.error.msg:
            self.__writer('%s\n' % self.error.msg.pop(0))
        self.error.col_offset = 0
        self.error.lineno = None
        self.error.raised = True


    def onError(self,exception,node,msg=None,with_tback=None):
        """ wrapper for raising exceptions from interpreter"""
        err_type,err_value,err_tback = sys.exc_info()
        
        # print("onError ",  node, msg) # self.dump(node,include_attributes=True))
        if self.error_msg is None: self.error_msg = []
        
        if self.error.node is None and node is not None:
            self.error.node = node
            # print('Errant node: ', self.dump(node))
            if isinstance(node, ast.Module):
                self.error.node = node.body[0]
            self.error.type = err_value.__class__.__name__
            self.error.msg.append("%s: %s" % (err_value.__class__.__name__, ' '.join(err_value.args)))

            if hasattr(node,'lineno'):     self.error.lineno = node.lineno
            if hasattr(node,'col_offset'): self.error.col_offset = node.col_offset
            
            if self.error.lineno is not None:
                self.lineno = self.lineno+self.error.lineno

            extra = ''
            if self.error.col_offset is not None:
                extra = ', column %i' % self.error.col_offset
            self.error.msg.append("%s, line %i%s" % (self.fname, self.lineno, extra))

            # print("Set error node to " ,self.error.node, self.error.lineno)

        if msg is not None:   self.error.msg.append(msg)

    
    def dump(self, node,**kw):  return ast.dump(node,**kw)

    # handlers for ast components
    def _NodeValue(self,node):  return self.interp(node.value)

    def doExpr(self,node):   return self._NodeValue(node)  # ('value',)
    def doIndex(self,node):  return self._NodeValue(node)  # ('value',)
    def doReturn(self,node): # ('value',)
        self.retval = self._NodeValue(node)
        return
    
    def doRepr(self,node):   return repr(self._NodeValue(node))  # ('value',)

    def doModule(self,node):    # ():('body',) 
        out = None
        for n in node.body: out = self.interp(n)
        return out

    def doExpression(self,node): return self.doModule(node) # ():('body',) 

    def doPass(self,node):    return None  # () 
    def doEllipsis(self,node): return Ellipsis # ??? 

    # for break and continue: set the instance variable _interrupt
    def doInterrupt(self,node):    # ()
        self._interrupt = node
        return node

    def doBreak(self,node):     return self.doInterrupt(node)
    def doContinue(self,node):  return self.doInterrupt(node)    

    def doAssert(self,node):    # ('test', 'msg')
        if not self.interp(node.test):
            self.onError(Assertiself.onError, node, msg=self.interp(node.msg()))
        return True

    def doList(self,node):    # ('elts', 'ctx') 
        return [self.interp(e) for e in node.elts]

    def doTuple(self,node):    # ('elts', 'ctx') 
        return tuple(self.doList(node))
    
    def doDict(self,node):    # ('keys', 'values')
        nodevals = zip(node.keys,node.values)
        interp = self.interp
        return dict([(interp(k),interp(v)) for k,v in nodevals])

    def doNum(self,node):  return node.n  # ('n',) 
    def doStr(self,node):  return node.s  # ('s',)

    def doName(self,node):    # ('id', 'ctx')
        """ Name node """
        ctx = _Context(node)
        if ctx == ast.Del:
            val = self.delSymbol(node.id)
        elif ctx == ast.Param:  # for Function Def
            val = str(node.id)
        else:
            val = self.getSymbol(node.id)
            if isinstance(val,DefinedVariable): val = val.evaluate()
        return val

    def _NodeAssign(self,n,val):
        """here we assign a value (not the node.value object) to a node
        this is used by doAssign, but also by for, list comprehension, etc.
        """
        if _Class(n) == ast.Name:
            sym = self.setSymbol(n.id,value=val)
            
        elif _Class(n) == ast.Attribute:
            if _Context(n) == ast.Load:
                self.onError(TDLError, n,
                        msg=  "canot Assign attribute %s???" % n.attr)
            setattr(self.interp(n.value),n.attr,val)

        elif _Class(n) == ast.Subscript:
            sym    = self.interp(n.value)
            slice  = self.interp(n.slice)
            if isinstance(n.slice,ast.Index):
                sym.__setitem__(slice,val)
            elif isinstance(n.slice,ast.Slice):
                sym.__setslice__(slice.start,slice.stop,val)
            elif isinstance(n.slice,ast.ExtSlice):
                sym[(slice)] = val
        elif _Class(n) in (ast.Tuple,ast.List):
            if len(val) == len(n.elts):
                for el,v in zip(n.elts,val):
                    self._NodeAssign(el,v)
            else:
                self.onError(ValueError, node, msg='too many values to unpack')

    def doAttribute(self,node):    # ('value', 'attr', 'ctx')
        ctx = _Context(node)
        # print("doAttribute",node.value,node.attr,ctx)
        if ctx == ast.Load:
            sym = self._NodeValue(node)
            if hasattr(sym,node.attr):
                val = getattr(sym,node.attr)
                if isinstance(val,DefinedVariable): val = val.evaluate()
                return val
            else:
                msg = "'%s' does not have an '%s' attribute" 
                self.onError(TDLError, node, msg= msg % (node.value,node.attr))

        elif ctx == ast.Del:
            return delattr(sym,attr)
        elif ctx == ast.Store:
            self.onError(TDLError,node,
                    msg = "attribute for storage: shouldn't be here!!")

    def doAssign(self,node):    # ('targets', 'value')
        val = self._NodeValue(node)
        for n in node.targets:  self._NodeAssign(n,val)
        return # return val

    def doAugAssign(self,node):    # ('target', 'op', 'value')
        return self.doAssign(ast.Assign(targets= [node.target],
                                        value  = ast.BinOp(left = node.target,
                                                           op   = node.op,
                                                           right= node.value)))
       
    def doSlice(self,node):    # ():('lower', 'upper', 'step')
        return slice(self.interp(node.lower),
                     self.interp(node.upper),
                     self.interp(node.step))

    def doExtSlice(self,node):    # ():('dims',)
        return tuple([self.interp(n) for n in node.dims])
    
    def doSubscript(self,node):    # ('value', 'slice', 'ctx') 
        # print("doSubscript: ", ast.dump(node))
        val    = self.interp(node.value)
        slice  = self.interp(node.slice)
        ctx = _Context(node)
        # print("doSubscript: ", val, slice, ctx)
        if ctx == ast.Load:
            if isinstance(node.slice,(ast.Index,ast.Slice)):
                return val.__getitem__(slice)
            elif isinstance(node.slice,ast.ExtSlice):
                return val[(slice)]
        elif ctx == ast.Store:
            self.onError(TDLError, node,
                    msg="subscript for storage: shouldn't be here!!")

    def doDelete(self,node):    # ('targets',)
        ctx = _Context(node)
        assert ctx == ast.Del, 'wrong Context for delete???'
        for n in node.targets: self.interp(n)

    def doUnaryOp(self,node):    # ('op', 'operand') 
        return _Op(node.op)(self.interp(node.operand))

    def doBinOp(self,node):    # ('left', 'op', 'right')
        return _Op(node.op)(self.interp(node.left), self.interp(node.right))

    def doBoolOp(self,node):    # ('op', 'values')
        val = self.interp(node.values.pop(0))
        for n in node.values:
            val =  _Op(node.op)(val,self.interp(n))
        return val
    
    def doCompare(self,node):    # ('left', 'ops', 'comparators')
        lval = self.interp(node.left)
        out  = True
        for op,rnode in zip(node.ops,node.comparators):
            rval = self.interp(rnode)
            out  = out and  _Op(op)(lval,rval)
            lval = rval
            if not out: break
        return out

    def doPrint(self,node):    # ('dest', 'values', 'nl') 
        """ note: implements Python2 style print statement, not
        print() function.  Probably, the 'tdl2py' translation
        should look for and translate print -> print_() to become
        a customized function call.
        """
        l = [self.interp(n) for n in  node.values]
        dest = self.interp(node.dest) or sys.stdout
        end = ''
        if node.nl: end = '\n'
        print(*l,file=dest,end=end)
        
    def doIf(self,node):    # ('test', 'body', 'orelse')
        block = node.orelse
        if self.interp(node.test): block = node.body
        for n in block: self.interp(n)

    def doIfExp(self,node):    # ('test', 'body', 'orelse') 
        expr = node.orelse
        if self.interp(node.test): expr = node.body
        return self.interp(exr)

    def doWhile(self,node):    # ('test', 'body', 'orelse')
        while self.interp(node.test):
            self._interrupt = None
            for n in node.body:
                v = self.interp(n)
                if self._interrupt is not None:  break
            if isinstance(self._interrupt,ast.Break): break
        else:
            for n in node.orelse: self.interp(n)
        self._interrupt = None

    def doFor(self,node):    # ('target', 'iter', 'body', 'orelse')
        for val in self.interp(node.iter):
            self._NodeAssign(node.target,val)
            self._interrupt = None
            for n in node.body:
                v = self.interp(n)
                if self._interrupt  is not None:  break
            if isinstance(self._interrupt,ast.Break): break
        else:
            for n in node.orelse: self.interp(n)
        self._interrupt = None

    def doListComp(self,node):    # ('elt', 'generators') 
        out = []
        for n in node.generators:
            if _Class(n) == ast.comprehension:
                for val in self.interp(n.iter):
                    self._NodeAssign(n.target,val)
                    add = True
                    for cond in n.ifs:
                        add = add and self.interp(cond)
                    if add:
                        out.append(self.interp(node.elt))
        return out

    def doCall(self,node):    # ('func', 'args', 'keywords', 'starargs', 'kwargs')
        func = self.interp(node.func)
        if not callable(func):
            self.onError(TDLError, node, msg="'%s' is not not callable" % (func))

        args = [self.interp(a) for a in node.args]
        if node.starargs is not None:
            args = args + self.interp(node.starargs)
        
        keywords = {}
        for k in node.keywords:
            if not isinstance(k,ast.keyword):
                self.onError(TDLError, node,
                        msg="keyword error in function call '%s'" % (func))
            keywords[k.arg] = self.interp(k.value)
        if node.kwargs is not None:  keywords.update(self.interp(node.kwargs))

        return func(*args,**keywords)
    
    def doFunctionDef(self,node):    # ('name', 'args', 'body', 'decorator_list') 
        if node.decorator_list != []:
            print("Warning: decorated procedures not supported!")
        args   = []
        kwargs = {}
        while node.args.defaults:
            defval = self.interp(node.args.defaults.pop())
            key    = self.interp(node.args.args.pop())
            kwargs[key] = defval
            
        for n in node.args.args: args.append(n.id)
     
        proc = Procedure(node.name, compiler=self,
                         body  =node.body,
                         args  =args,
                         kwargs=kwargs,
                         vararg=node.args.vararg,
                         varkws=node.args.kwarg)
        self.setSymbol(node.name,value=proc)

    def doImport(self,node):    # ('names',) 
        for n in node.names:
            self.import_module(n.name,asname=n.asname)
        
    def doImportFrom(self,node):    # ('module', 'names', 'level') 
        fromlist, asname = [], []
        for n in node.names:
            fromlist.append(n.name)
            asname.append(n.asname)
        self.import_module(node.module,asname=asname,fromlist=fromlist)

    ## 
    # not yet implemented:
    ## 
    def doExceptHandler(self,node): # ('type', 'name', 'body')
        print("except handler")
        return None
    
    def doTryExcept(self,node):    # ('body', 'handlers', 'orelse') 
        print("Incomplete Try Except")
        print( self.dump(node))
        for n in node.body:
            try:
                self.interp(n)
            except:
                e_type,e_value,e_tback = sys.exc_info()
                handled = False
                for h in node.handlers:
                    handler_type = self.interp(h.type)
                    if handler_type is None or isinstance(e_type(),handler_type):
                        handled=True
                        self._NodeAssign(h.name,e_value)
                        for b in h.body: self.interp(b)
                        break
                if not handled:
                    # print("unhandled exception: ", dir(e_tback), e_tback.tb_lineno, e_tback.tb_lasti)
                    self.onError(e_type, n, msg=e_value) # print("%s: %s" % (e_type.__name__, e_value))
                    
    def doTryFinally(self,node):    # ('body', 'finalbody') 
        return self.NotImplemented(node)
        
    def doGeneratorExp(self,node):    # ('elt', 'generators') 
        print('Incomplete GeneratorExp ')
        print(ast.dump(node.elt),include_attributes=True)
        for n in node.generators:
            print(n)             

    def doRaise(self,node):    # ('type', 'inst', 'tback') 
        print(self.dump(node,include_attributes=True))
        self.onError(TDLError, node)
        #, self.interp(node.type),self.interp(node.inst),tback=self.interp(node.tback))



    def import_module(self, name, asname=None, fromlist=None, reload=False):
        """
        import a module (tdl or python), installing it into the symbol table.
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

        this method covers a lot of cases (tdl or python, import
        or from-import, use of asname) and so is fairly long.
        """
        symtable = self.symtable
        _sys     = symtable._sys

        # step 1  import the module to a global location
        #   either sys.modules for python modules
        #   or  _sys.modules for tdl modules
        # reload takes effect here in the normal python way:
        #   if
        if reload or (name not in sys.modules and
                      name not in _sys.modules):
            
            # first look for "name.tdl"
            isTDL = False
            tdlname = "%s.tdl" % name
            for dirname in _sys.path:
                if tdlname in os.listdir(dirname):
                    isTDL = True
                    modname = os.path.abspath(os.path.join(dirname,tdlname))

                    # save current module group
                    #  create new group, set as moduleGroup and localGroup
                    saveGroups = _sys.localGroup,_sys.moduleGroup
                    thismod = symbolTable.Group()
                    _sys.modules[name]= thismod
                    symtable._set_local_mod((thismod,thismod))

                    text = open(modname).read()
                    inptext = inputText.InputText()
                    inptext.put(text,filename=modname)
                       
                    while inptext:
                        block,fname,lineno = inptext.get()
                        self.eval(block,fname=fname,lineno=lineno)

                        while self.error.msg:
                            self.__writer('%s\n' % self.error.msg.pop(0))


                        if self.error.raised:
                            print("..break", len(inptext))
                            break
                        
                    thismod = _sys.modules[name]               
                    symtable._set_local_mod(saveGroups)
            if self.error_raised:
                print("import returning")
                return
                        

            # or, if not a tdl module, load as a regular python module
            if not isTDL:
                __import__(name)
                thismod = sys.modules[name]
               
        else: # previously loaded module, just do lookup
            if name in _sys.modules:
                thismod = _sys.modules[name]
            elif name in sys.modules:
                thismod = sys.modules[name]               
               
        # now we install thismodule into the current moduleGroup
        # import full module
        if fromlist is None:
            if asname is None: asname = name
            parts = asname.split('.')
            asname = parts.pop()
            top = _sys.moduleGroup
            while len(parts)>0:
                subname = parts.pop(0)
                subgrp  = Group()
                setattr(top, subname, subgrp)
                top = subgrp
                
            setattr(top, asname, thismod)
        # import-from construct
        else:
            if asname is None:  asname = [None]*len(fromlist)
            for sym,alias in zip(fromlist,asname):
                if alias is None: alias = sym
                setattr(_sys.moduleGroup, alias, getattr(thismod,sym))
    # end
