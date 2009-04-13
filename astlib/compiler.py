from __future__ import division, print_function
import os
import sys
import ast
from   symbolTable import symbolTable
import __builtin__
import numpy
import copy

class TDLError(Exception):
    def __init__(self,error,descr = None,node = None):
        self.error = error
        self.descr = descr
    def __repr__(self):
        return "%s: %s" % (self.error, self.descr)
    __str__ = __repr__

def onError(exception,msg,tback=None):
    """ wrapper for raising exceptions from interpreter"""
    err_type,err_value,err_tback = sys.exc_info()
    print(" onError:  exc_info  =  ", sys.exc_info())
    print(" onError:  exception =  ", exception)
    print(" onError:  message   =  ", msg)
    
    raise exception, msg

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
def _Op(node):        return _operators[_Class(node)]


# inherit these from python's __builtin__
_from_builtins= ('ArithmeticError', 'AssertionError', 'AttributeError',
                 'BaseException', 'BufferError', 'BytesWarning',
                 'DeprecationWarning', 'EOFError', 'EnvironmentError',
                 'Exception', 'False', 'FloatingPointError',
                 'GeneratorExit', 'IOError', 'ImportError', 'ImportWarning',
                 'IndentationError', 'IndexError', 'KeyError',
                 'KeyboardInterrupt', 'LookupError', 'MemoryError',
                 'NameError', 'None', 'NotImplemented',
                 'NotImplementedError', 'OSError', 'OverflowError',
                 'ReferenceError', 'RuntimeError', 'RuntimeWarning',
                 'StandardError', 'StopIteration', 'SyntaxError',
                 'SyntaxWarning', 'SystemError', 'SystemExit', 'True',
                 'TypeError', 'UnboundLocalError', 'UnicodeDecodeError',
                 'UnicodeEncodeError', 'UnicodeError',
                 'UnicodeTranslateError', 'UnicodeWarning', 'ValueError',
                 'Warning', 'ZeroDivisionError', 'abs', 'all', 'any',
                 'apply', 'basestring', 'bin', 'bool', 'buffer',
                 'bytearray', 'bytes', 'callable', 'chr', 'cmp', 'coerce',
                 'complex', 'delattr', 'dict', 'dir', 'divmod', 'enumerate',
                 'file', 'filter', 'float', 'format', 'frozenset',
                 'getattr', 'hasattr', 'hash', 'hex', 'id', 'int',
                 'isinstance', 'len', 'list', 'map', 'max', 'min', 
                 'oct', 'open', 'ord', 'pow', 'property', 'range',
                 'raw_input', 'reduce', 'repr', 'reversed', 'round', 'set',
                 'setattr', 'slice', 'sorted', 'str', 'sum', 'tuple',
                 'type', 'unichr', 'unicode', 'zip')

# inherit these from numpy
_from_numpy = ('pi','e', 'array','sin','cos','tan','exp','log','log10',
               'sqrt','arange', 'arccos', 'arccosh', 'arcsin', 'arcsinh',
               'arctan', 'arctan2', 'arctanh', 'argmax', 'argmin',
               'argsort', 'array', 'cosh', 'fabs', 'floor', 'floor_divide',
               'fmod', 'tanh', 'sign', 'sinh', 'identity', 'take',
               'choose', 'add', 'allclose', 'alltrue', 'around', 'asarray',
               'average', 'bitwise_and', 'bitwise_or', 'bitwise_xor',
               'ceil', 'clip', 'compress', 'concatenate', 'conjugate',
               'convolve', 'cumproduct', 'cumsum', 'diagonal', 'divide',
               'dot', 'equal', 'greater', 'greater_equal', 'hypot',
               'indices', 'invert', 'left_shift', 'less', 'less_equal',
               'logical_and', 'logical_not', 'logical_or', 'logical_xor',
               'maximum', 'minimum', 'multiply', 'negative', 'nonzero',
               'not_equal', 'ones', 'power', 'product', 'put', 'putmask',
               'rank', 'ravel', 'remainder', 'repeat', 'reshape', 'resize',
               'right_shift', 'searchsorted', 'shape', 'size', 'sometrue',
               'sort', 'subtract', 'sum', 'swapaxes', 'trace', 'transpose',
               'true_divide', 'vdot', 'where', 'zeros')

def group(compiler=None,**kw):
    try:
        g = compiler.symtable.createGroup()
        for k,v in kw.items():  setattr(g,k,v)
        return g
    except:
        return None

def showgroup(gname,compiler=None):
    if compiler is not None:
        compiler.symtable.show_group(gname)

def definevar(name,expr,compiler=None):
    # print("===DEFVAR ", name, expr, compiler)
    if compiler is not None:
        defvar = DefinedVariable(expr=expr,compiler=compiler)
        compiler.symtable.setSymbol(name,defvar)


def _copy(obj,**kw):
    if kw.has_key('compiler'):
        compiler = kw.pop('compiler')
    return copy.deepcopy(obj)

        
_local_funcs = {'group':group,
                'showgroup':showgroup,
                'definevar':definevar,
                'copy': _copy,
                }
        
####
class DefinedVariable(object):
    """defined variable: re-evaluate on access

    Note that the LocalGroup/ModuleGroup are cached
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
            self._groups = (_sys.LocalGroup,_sys.ModuleGroup)

    def evaluate(self):
        if self.ast is None: self.compile()

        if self.ast is None:
            msg = "Cannot compile '%s'"  % (self.expr)
            raise Warning, msg            
            
        if hasattr(self.compiler,'interp'):
            _sys = self.compiler.symtable._sys
            # save current LocalGroup/ModuleGroup 
            save_groups  = _sys.LocalGroup,_sys.ModuleGroup
            
            _sys.LocalGroup,_sys.ModuleGroup = self._groups
            rval = self.compiler.interp(self.ast)

            _sys.LocalGroup,_sys.ModuleGroup = save_groups
            return rval
        else:
            msg = "Cannot evaluate '%s'"  % (self.expr)
            raise ValueError, msg

     
class Procedure(object):
    """tdl procedure:  function """
    def __init__(self, name, compiler=None,
                 body=None, args=None, kwargs=None,
                 vararg=None, varkws=None):
        self.name     = name
        self.compiler = compiler
        self.module   = compiler.symtable._sys.ModuleGroup
        self.body     = body
        self.argnames = args
        self.kwargs   = kwargs
        self.vararg   = vararg
        self.varkws   = varkws
        
    def __call__(self,*args,**kwargs):
        stable  = self.compiler.symmtable
        sys     = stable._sys      
        lgroup  = stable.createGroup()
        locname = '%s_%s' % (self.name,hex(id(lgroup))[2:])

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
        
        setattr(stable,locname,lgroup)
                
        grps_save = sys.LocalGroup,sys.ModuleGroup
        stable._set_local_mod((locname, self.module))
        
        retval = None
        self.compiler.retval = None
        for node in self.body:
            self.compiler.interp(node)
            if self.compiler.retval is not None:
                retval = self.compiler.retval
                break

        delattr(stable,locname)
        sys.LocalGroup,sys.ModuleGroup = grps_save
        stable._set_local_mod(grps_save)
        self.compiler.retval = None
        return retval
    
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
    def __init__(self,symtable=None,input=None, output=None,libs=None):
        # if symtable is None:
        symtable = symbolTable()
        self.setSymbol   = symtable.setSymbol
        self.getSymbol   = symtable.getSymbol
        self.delSymbol   = symtable.delSymbol        
        self._interrupt  = None
        self.symtable    = symtable
        self.retval      = None

        self.setSymbol('_builtin.None',None)
        load = symtable._load_functions
        load(_from_builtins, group=symtable._builtin, parent=__builtin__)        
        load(_from_numpy,  group=symtable._math, parent=numpy)
        load(_local_funcs, group=symtable._builtin, compiler=self)

    ##
    def NotImplemented(self,node):
        onError(TDLError, "'%s' not supported" % (_ClassName(node)))

    # main entry point for Ast node evaluation
    #  compile:  string statement -> ast
    #  interp :  ast -> result
    #  eval   :  string statement -> result = interp(compile(statement))
    def compile(self,text):
        """compile statement/expression to Ast representation    """
        try:
            return ast.parse(text)
        except:
            onError(SyntaxError, text)
            
    
    def dump(self, node,**kw):  return ast.dump(node,**kw)
        
    def interp(self, node):
        """executes compiled Ast representation for an expression"""
        # it **is** important to keep this, as
        # internal code here may run interp(None)
        # and expect a None in return.
        if node is None: return None
        if isinstance(node,(str,unicode)):  node = self.compile(node)

        methodName = "do%s" % _Class(node).__name__
        return getattr(self,methodName)(node)

    def eval(self,expr):
        """evaluates a single statement"""
        return self.interp(self.compile(expr))


    # handlers for ast components
    def doModule(self,node):    # ():('body',) 
        out = None
        for n in node.body: out = self.interp(n)
        return out

    def doExpression(self,node): return self.doModule(node) # ():('body',) 

    def _NodeValue(self,node):
        'common case'
        return self.interp(node.value)

    def doExpr(self,node):   return self._NodeValue(node)  # ('value',)
    def doIndex(self,node):  return self._NodeValue(node)  # ('value',)
    def doReturn(self,node): # ('value',)
        self.retval = self._NodeValue(node)
        return
    
    def doRepr(self,node):   return repr(self._NodeValue(node))  # ('value',)

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
            onError(AssertionError, self.interp(node.msg()))
        return True

    def doList(self,node):    # ('elts', 'ctx') 
        return [self.interp(e) for e in node.elts]

    def doTuple(self,node):    # ('elts', 'ctx') 
        return tuple(self.doList(node))
    
    def doDict(self,node):    # ('keys', 'values') 
        return dict([(self.interp(k),self.interp(v)) for k,v in zip(node.keys,node.values)] )

    def doNum(self,node):  return node.n  # ('n',) 
    def doStr(self,node):  return node.s  # ('s',)

    

    def doName(self,node):    # ('id', 'ctx')
        """ Name node """
        ctx = _Context(node)
        method = self.getSymbol
        if ctx == ast.Del: method = self.delSymbol
        if ctx == ast.Param: method = str  # for Function Def
        val = method(node.id)
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
                onError(TDLError,
                        "canot Assign attribute %s???" % n.attr)
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
                onError(ValueError, 'too many values to unpack')

    def doAttribute(self,node):    # ('value', 'attr', 'ctx')
        ctx = _Context(node)
        # print("doAttribute",node.value,node.attr,ctx)
        if ctx == ast.Load:
            sym = self._NodeValue(node)
            if hasattr(sym,node.attr):
                return getattr(sym,node.attr)
            else:
                msg = "'%s' does not have an '%s' attribute" 
                onError(TDLError, msg % (node.value,node.attr))

        elif ctx == ast.Del:
            return delattr(sym,attr)
        elif ctx == ast.Store:
            onError(TDLError,
                    "attribute for storage: shouldn't be here!!")

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
            onError(TDLError,
                    "subscript for storage: shouldn't be here!!")

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

    def doCall(self,node):    # ('func', 'args', 'keywords', 'starargs', 'kwargs')
        func = self.interp(node.func)
        if not callable(func):
            onError(TDLError, "'%s' is not not callable" % (func))

        args = [self.interp(a) for a in node.args]
        if node.starargs is not None:
            args = args + self.interp(node.starargs)
        
        keywords = {}
        for k in node.keywords:
            if not isinstance(k,ast.keyword):
                onError(TDLError,
                        "keyword error in function call '%s'" % (func))
            keywords[k.arg] = self.interp(k.value)
        if node.kwargs is not None:  keywords.update(self.interp(node.kwargs))

        return func(*args,**keywords)
    
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
                    onError(e_type, e_value,e_tback) # print("%s: %s" % (e_type.__name__, e_value))
                    
    def doTryFinally(self,node):    # ('body', 'finalbody') 
        return self.NotImplemented(node)
        
    def doGeneratorExp(self,node):    # ('elt', 'generators') 
        print('Incomplete GeneratorExp ')
        print(ast.dump(node.elt),include_attributes=True)
        for n in node.generators:
            print(n)             

    def doRaise(self,node):    # ('type', 'inst', 'tback') 
        print(" Hello from Raise!" )
        print(self.dump(node,include_attributes=True))
        onError(self.interp(node.type),self.interp(node.inst),tback=self.interp(node.tback))

    def doImport(self,node):    # ('names',) 
        for n in node.names:
            print(self.dump(n))
            name,asname = n.name, n.asname
            self.symtable.import_module(name,asname=asname)
        # return self.NotImplemented(node)
        
    def doImportFrom(self,node):    # ('module', 'names', 'level') 
        print(":: ", node.module, node.level)
        self.symtable.import_module(node.module)
        for n in node.names:
            print("Dump: ", self.dump(n))
            name,asname = n.name, n.asname
            # self.symtable.import_module(name,asname=asname)


    #
