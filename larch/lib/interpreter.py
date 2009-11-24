from __future__ import division, print_function
import os
import sys
import ast
import traceback
import inputText
import symbolTable
import builtins
import numpy
import __builtin__
from util import closure

__version__ = '0.9.1'
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

##
class DefinedVariable(object):
    """defined variable: re-evaluate on access

    Note that the localGroup/moduleGroup are cached
    at compile time, and restored for evaluation.
    """
    def __init__(self, expr=None, larch=None):
        self.expr = expr
        self.larch = larch
        self.ast = None
        self._groups = None,None
        self.compile()

    def __repr__(self):
        return "<DefinedVariable: '%s'>" % (self.expr)
        
    def compile(self):
        if self.larch is not None and self.expr is not None:
            self.ast = self.larch.compile(self.expr)
            _sys = self.larch.symtable._sys
            self._groups = (_sys.localGroup,_sys.moduleGroup)

    def evaluate(self):
        if self.ast is None: self.compile()

        if self.ast is None:
            msg = "Cannot compile '%s'"  % (self.expr)
            raise Warning(msg)
            
        if hasattr(self.larch,'interp'):
            _sys = self.larch.symtable._sys
            # save current localGroup/moduleGroup 
            save_groups  = _sys.localGroup,_sys.moduleGroup
            
            _sys.localGroup,_sys.moduleGroup = self._groups
            rval = self.larch.interp(self.ast,expr=self.expr)

            _sys.localGroup,_sys.moduleGroup = save_groups
            return rval
        else:
            msg = "Cannot evaluate '%s'"  % (self.expr)
            raise ValueError, msg

class Procedure(object):
    """larch procedure:  function """
    def __init__(self, name, larch=None, doc=None,
                 fname='<StdIn>', lineno=0,
                 body=None, args=None, kwargs=None,
                 vararg=None, varkws=None):
        self.name     = name
        self.larch = larch
        self.modgroup = larch.symtable._sys.moduleGroup
        self.body     = body
        self.argnames = args
        self.kwargs   = kwargs
        self.vararg   = vararg
        self.varkws   = varkws
        self.__doc__  = doc
        self.lineno   = lineno
        self.fname    = fname
        
    def __repr__(self):
        sig = "%s(" % self.name
        if len(self.argnames )>0:
            sig = "%s%s" % (sig,','.join(self.argnames))
        if len(self.kwargs)>0:
            if len(self.argnames )>0:
                sig = "%s," % sig            
            _kw = ["%s=%s" % (k,v) for k,v in self.kwargs.items()]
            sig = "%s%s" % (sig,','.join(_kw))
            
        sig = "<Procedure %s, file=%s)>" % (sig,self.fname)

        if self.__doc__ is not None:
            sig = "%s\n  %s" % (sig, self.__doc__)
        return sig

    def __call__(self,*args,**kwargs):
        try:
            stable  = self.larch.symtable
            _sys    = stable._sys      
            lgroup  = stable.createGroup()
        except:
            self.larch.addException(None,
                                    msg='Cannot run Procedure %s'%self.name,
                                    expr='<>',
                                    fname=self.fname,lineno=self.lineno+1,
                                    py_exc=sys.exc_info())

        args   = list(args)
        n_args = len(args)
        n_expected = len(self.argnames)
        if n_args < n_expected:
            msg='not enough arguments for Procedure %s (expected %i, got %i)'%(self.name,
                                                                               n_expected,
                                                                               n_args)
            self.larch.addException(None,msg=msg, expr='<>',
                                    fname=self.fname,lineno=self.lineno+1,
                                    py_exc=sys.exc_info())
                
        for argname in self.argnames:
            setattr(lgroup, argname,args.pop(0))
            
        try:         
            if self.vararg is not None and len(args)>0:
                setattr(lgroup, self.vararg, args)
                
            for k,v in self.kwargs.items():
                if kwargs.has_key(k):  v = kwargs.pop(k)
                setattr(lgroup, k, v)

            if self.varkws is not None:
                setattr(lgroup, self.varkws, kwargs)

            if len(kwargs) > 0:
                msg='extra keyword arguments for Procedure %s (%s)'%(self.name,
                                                                     ','.join(kwargs.keys()))
                self.larch.addException(None,msg=msg, expr='<>',
                                        fname=self.fname,lineno=self.lineno+1,
                                        py_exc=sys.exc_info())

                
        except:
            self.larch.addException(None,
                                    msg='incorrect arguments for Procedure %s'%self.name,
                                    expr='<>',
                                    fname=self.fname,lineno=self.lineno+1,
                                    py_exc=sys.exc_info())            
            

        grps_save = _sys.localGroup,_sys.moduleGroup
        stable._set_local_mod((lgroup, self.modgroup))
        
        retval = None
        self.larch.retval = None
        # print("***** Calling proc " , self.name, self.fname, self.lineno)
        for node in self.body:
            self.larch.interp(node,expr='<>',fname=self.fname,lineno=self.lineno)
            
            if len(self.larch.error)>0:
                break
            if self.larch.retval is not None:
                retval = self.larch.retval
                break

        _sys.localGroup,_sys.moduleGroup = grps_save
        stable._set_local_mod(grps_save)
        self.larch.retval = None
        del lgroup
        return retval
    
class LarchExceptionHolder:
    def __init__(self,node,msg='',fname='<StdIn>',
                 py_exc=(None,None),
                 expr=None, lineno=0):
        self.node   = node
        self.fname  = fname
        self.expr   = expr
        self.msg    = msg
        self.py_exc = py_exc
        self.lineno = lineno
        self.exc_info = sys.exc_info()

    def get_error(self):
        node = self.node
        node_lineno = 0
        node_col_offset = 0

        if node is not None:
            try:
                node_lineno = node.lineno
                node_col_offset = self.node.col_offset
            except: 
                pass
            
        lineno = self.lineno + node_lineno 

        exc_text = str(self.exc_info[1])
        if exc_text in (None,'None'): exc_text = ''
        expr = self.expr
        if expr == '<>': # denotes non-saved expression -- go fetch from file!
            ftmp = open(self.fname,'r')
            expr = ftmp.readlines()[lineno-1][:-1]
            ftmp.close()

        out = []
        if len(exc_text) > 0: out.append(exc_text)
        py_etype, py_eval = self.py_exc
        
        #         if py_etype is not None and py_eval is not None:
        #             out.append("%s: %s" % (py_etype, py_eval))
        if (self.fname == '<StdInput>' and self.lineno==0):
            out.append(' ')
        else:
            out.append(" %s, line number %i" % (self.fname,self.lineno))
            
        out.append("     %s" % expr)
        if node_col_offset>0:
            out.append("    %s^^^" % ((node_col_offset)*' '))

        return (self.msg,'\n'.join(out))


class LarchError(Exception):
    def __init__(self,error='Larch Error', descr=None, node=None):
        self.error = error
        self.descr = descr
        self.node  = node
    def __repr__(self): return "%s: %s" % (self.error, self.descr)
    __str__ = __repr__




class Interpreter:
    """larch program compiler and interpreter.
  This module compiles expressions and statements to AST representation,
  using python's ast module, and then executes the AST representation
  using a custom SymboplTable for named object (variable, functions).
  This then gives a restricted version of Python, with slightly modified
  namespace rules.  The program syntax here is expected to be valid Python,
  but that may have been translated as with the inputText module.

    
  The following Python syntax is not supported:
      Exec, Lambda, Class, Global, Generators, Yield, Decorators
        
  In addition, Function is greatly altered so as to allow a Larch procedure.
    
    """
    def __init__(self,symtable=None,input=None, writer=None):

        self.__writer = writer or sys.stdout.write
       
        if symtable is None:
            symtable = symbolTable.symbolTable()
        self.isGroup    = symtable.isGroup
        self.symtable   = symtable
        self.setSymbol  = symtable.setSymbol
        self.getSymbol  = symtable.getSymbol
        self.delSymbol  = symtable.delSymbol        
        self._interrupt = None
        self.error      = []
        self.expr       = None
        self.retval     = None

        imports = ((builtins._from_builtin, __builtin__,'_builtin'),
                   (builtins._from_numpy,   numpy, '_math'))
        
        for symlist, pymod, larchmod in imports:
            group = getattr(symtable,larchmod)
            for sym in symlist:
                setattr(group, sym, getattr(pymod, sym))

        group = getattr(symtable, '_builtin')
        for fname,fcn in builtins._local_funcs.items():
            setattr(group, fname, closure(func=fcn,larch=self))
        setattr(group, 'definevar', closure(func=self.__definevar))

        
    def __definevar(self,name,expr):
        """define a defined variable (re-evaluate on access)"""
        defvar = DefinedVariable(expr=expr,larch=self)
        self.setSymbol(name,defvar)

    def NotImplemented(self,node):
        cname = node.__class__.__name__
        self.addException(node,"'%s' not supported" % (cname),
                          py_exc=sys.exc_info())
    
    
    # main entry point for Ast node evaluation
    #  compile:  string statement -> ast
    #  interp :  ast -> result
    #  eval   :  string statement -> result = interp(compile(statement))
                      
    def addException(self,node,msg='',expr=None,fname=None,lineno=0,py_exc=None):

        if self.error is None: self.error = []
        if expr  is None: expr  = self.expr
        if fname is None: fname = self.fname        
        if lineno is None: lineno = self.lineno

        if len(self.error) > 0 and not isinstance(node, ast.Module):
            msg = 'Extra Error (%s)' % msg

        if py_exc is None:
            etype,evalue = None,None
        else:
            etype,evalue,tb = py_exc

        # print(" ADD Except:", node, msg, expr, fname, lineno)
        err = LarchExceptionHolder(node,
                                   py_exc=(etype,evalue),
                                   msg=msg,
                                   expr= expr,
                                   fname= fname,
                                   lineno=lineno)
        self._interrupt = ast.Break()
           
            
        self.error.append(err)
        
    def compile(self,text,fname=None,lineno=0):
        """compile statement/expression to Ast representation    """
        self.expr  = text
        # print(" larch compile: '%s'" % text)
        try:
            return  ast.parse(text)

        except:
            self.addException(None,msg='Syntax Error',
                              expr=text,
                              fname=fname,
                              lineno=lineno,
                              py_exc=sys.exc_info())
            
            
        
    def interp(self, node, expr=None, fname=None, lineno=None):
        """executes compiled Ast representation for an expression"""

        # Note: keep the 'node is None' test: internal code here may run
        #    interp(None) and expect a None in return.
        if node is None: return None
        if isinstance(node,(str,unicode)):  node = self.compile(node)

               
        method = "do%s" % node.__class__.__name__
        if lineno is not None: self.lineno = lineno
        if fname  is not None: self.fname  = fname
        if expr   is not None: self.expr   = expr
        
        # print(" Interp: ", node, node.__class__.__name__,  fname, lineno)
        # print(" Interp: ", expr, node, method)
        ret = None
        try:
            fcn = getattr(self,method)
        except:
            self.addException(node,msg='Lookup Error',
                              expr=expr,fname=fname,lineno=lineno,
                              py_exc=sys.exc_info())

            return ret
        
        try:
            ret = fcn(node)
        except:
            self.addException(node,msg='Runtime Error',
                              expr=expr,fname=fname,lineno=lineno,
                              py_exc=sys.exc_info())              
            
        if isinstance(ret,enumerate):  ret = list(ret)
        return ret
            
    def eval(self,expr,fname=None,lineno=0):
        """evaluates a single statement"""
        self.fname = fname        
        self.lineno = lineno
        self.error = []
        node = self.compile(expr,fname=fname,lineno=lineno)
        if not self.error:
            o = self.interp(node,expr=expr,fname=fname,lineno=lineno)
            if len(self.error) > 0:
                self.addException(node,msg='Eval Error',
                                  expr=expr,fname=fname,lineno=lineno,
                                  py_exc=sys.exc_info())
                
                return o
            return o
            
        
    def dump(self, node,**kw):  return ast.dump(node,**kw)

    # handlers for ast components
    def doExpr(self,node):   return self.interp(node.value)  # ('value',)
    def doIndex(self,node):  return self.interp(node.value)  # ('value',)
    def doReturn(self,node): # ('value',)
        self.retval = self.interp(node.value)
        return
    
    def doRepr(self,node):   return repr(self.interp(node.value))  # ('value',)

    def doModule(self,node):    # ():('body',) 
        out = None
        for n in node.body: out = self.interp(n)
        return out

    def doExpression(self,node): return self.doModule(node) # ():('body',) 

    def doPass(self,node):    return None  # () 
    def doEllipsis(self,node): return Ellipsis

    # for break and continue: set the instance variable _interrupt
    def doInterrupt(self,node):    # ()
        self._interrupt = node
        return node

    def doBreak(self,node):     return self.doInterrupt(node)
    def doContinue(self,node):  return self.doInterrupt(node)    

    def doAssert(self,node):    # ('test', 'msg')
        if not self.interp(node.test):
            raise AssertError( self.interp(node.msg()) )
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
        ctx = node.ctx.__class__
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
        if n.__class__ == ast.Name:
            sym = self.setSymbol(n.id,value=val)
            
        elif n.__class__ == ast.Attribute:
            if n.ctx.__class__  == ast.Load:
                self.addException(node,
                                  "cannot assign to attribute %s" % n.attr)

            setattr(self.interp(n.value),n.attr,val)

        elif n.__class__ == ast.Subscript:
            sym    = self.interp(n.value)
            slice  = self.interp(n.slice)
            if isinstance(n.slice,ast.Index):
                sym.__setitem__(slice,val)
            elif isinstance(n.slice,ast.Slice):
                sym.__setslice__(slice.start,slice.stop,val)
            elif isinstance(n.slice,ast.ExtSlice):
                sym[(slice)] = val
        elif n.__class__ in (ast.Tuple,ast.List):
            if len(val) == len(n.elts):
                for el,v in zip(n.elts,val):
                    self._NodeAssign(el,v)
            else:
                raise ValueError('too many values to unpack')

    def doAttribute(self,node):    # ('value', 'attr', 'ctx')
        ctx = node.ctx.__class__
        # print("doAttribute",node.value,node.attr,ctx)
        if ctx == ast.Load:
            sym = self.interp(node.value)
            if hasattr(sym,node.attr):
                val = getattr(sym,node.attr)
                if isinstance(val,DefinedVariable): val = val.evaluate()
                return val
            else:
                obj = self.interp(node.value)
                fmt = "%s does not have member '%s'"                
                if not self.isGroup(obj):
                    obj = obj.__class__
                    fmt = "%s does not have attribute '%s'"

                msg= fmt % (obj, node.attr)                    
                self.addException(node,msg=msg, py_exc=sys.exc_info())

        elif ctx == ast.Del:
            return delattr(sym,attr)
        elif ctx == ast.Store:
            msg = "attribute for storage: shouldn't be here!"
            self.addException(node,msg=msg, py_exc=sys.exc_info())        

    def doAssign(self,node):    # ('targets', 'value')
        val = self.interp(node.value)
        for n in node.targets:  self._NodeAssign(n,val)
        return # return val

    def doAugAssign(self,node):    # ('target', 'op', 'value')
        # print( "AugASSIGN ", node.target, node.value)
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
        nslice = self.interp(node.slice)
        ctx = node.ctx.__class__
        if ctx in ( ast.Load, ast.Store):
            if isinstance(node.slice,(ast.Index,ast.Slice,ast.Ellipsis)):
                return val.__getitem__(nslice)
            elif isinstance(node.slice,ast.ExtSlice):
                return val[(nslice)]
        else:
            msg = "subscript with unknown context"
            self.addException(node,msg=msg, py_exc=sys.exc_info())
            

    def doDelete(self,node):    # ('targets',)
        for n in node.targets:
            ctx = n.ctx.__class__
            assert ctx == ast.Del, 'wrong Context for delete???'

            print(ast.dump(n))
            print(" READY TO DELETE:: ", self.interp(n))

    def doUnaryOp(self,node):    # ('op', 'operand') 
        op = _operators[node.op.__class__]
        return op(self.interp(node.operand))
    
    def doBinOp(self,node):    # ('left', 'op', 'right')
        op = _operators[node.op.__class__]        
        return op(self.interp(node.left), self.interp(node.right))

    def doBoolOp(self,node):    # ('op', 'values')
        val = self.interp(node.values.pop(0))
        op = _operators[node.op.__class__]
        isAnd = True
        if ast.Or == node.op.__class__: isAnd = False
        if (isAnd and val) or (not isAnd and not val):
            for n in node.values:
                val =  op(val,self.interp(n))
                if isAnd:
                    if (not val):  break
                else:
                    if val: break
        return val
    
    def doCompare(self,node):    # ('left', 'ops', 'comparators')
        lval = self.interp(node.left)
        out  = True
        for op,rnode in zip(node.ops,node.comparators):
            comp = _operators[op.__class__]            
            rval = self.interp(rnode)
            out  = out and  comp(lval,rval)
            lval = rval
            if not out: break
        return out

    def doPrint(self,node):    # ('dest', 'values', 'nl') 
        """ note: implements Python2 style print statement, not
        print() function.  Probably, the 'larch2py' translation
        should look for and translate print -> print_() to become
        a customized function call.
        """
        dest = self.interp(node.dest) or sys.stdout
        end = ''
        if node.nl: end = '\n'
        out = [self.interp(n) for n in node.values]
        if out and len(self.error)==0:
            print(*out, file=dest, end=end)
        
    def doIf(self,node):    # ('test', 'body', 'orelse')
        block = node.orelse
        if self.interp(node.test): block = node.body
        for n in block: self.interp(n)

    def doIfExp(self,node):    # ('test', 'body', 'orelse') 
        expr = node.orelse
        if self.interp(node.test): expr = node.body
        return self.interp(expr)

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
                if self._interrupt is not None:  break
            if isinstance(self._interrupt,ast.Break): break
        else:
            for n in node.orelse: self.interp(n)
        self._interrupt = None

    def doListComp(self,node):    # ('elt', 'generators') 
        out = []
        for n in node.generators:
            if n.__class__ == ast.comprehension:
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
            msg = "'%s' is not not callable" % (func)
            self.addException(node,msg=msg, py_exc=sys.exc_info())

        args = [self.interp(a) for a in node.args]
        if node.starargs is not None:
            args = args + self.interp(node.starargs)
        
        keywords = {}
        for k in node.keywords:
            if not isinstance(k,ast.keyword):
                msg = "keyword error in function call '%s'" % (func)
                self.addException(node,msg=msg, py_exc=sys.exc_info())
            
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
     
        doc = None
        if isinstance(node.body[0],ast.Expr):
            docnode = node.body.pop(0)
            doc = self.interp(docnode.value)

        proc = Procedure(node.name, larch=self,
                         doc=doc,
                         lineno = self.lineno,
                         body = node.body,
                         fname = self.fname,
                         args = args,
                         kwargs = kwargs,
                         vararg = node.args.vararg,
                         varkws = node.args.kwarg)
        self.setSymbol(node.name,value=proc)

    # imports
    def doImport(self,node):    # ('names',) 
        for n in node.names:
            self.import_module(n.name,asname=n.asname)
        
    def doImportFrom(self,node):    # ('module', 'names', 'level') 
        fromlist, asname = [], []
        for n in node.names:
            fromlist.append(n.name)
            asname.append(n.asname)
        self.import_module(node.module,asname=asname,fromlist=fromlist)


    def import_module(self, name, asname=None, fromlist=None, reload=False):
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
        _sys     = symtable._sys
        
        # print("  larch path: ", _sys.path)

        for i in _sys.path:
            if i not in sys.path:
                sys.path.append(i)

        # print("  python path: ", sys.path)
        
        # step 1  import the module to a global location
        #   either sys.modules for python modules
        #   or  _sys.modules for larch modules
        # reload takes effect here in the normal python way:
        #   if
        do_load = (name not in _sys.modules and name not in sys.modules)
        do_load = do_load or reload
        

        # print("  larch do_load ", do_load)
        if do_load:
            # first look for "name.lar"
            isLarch = False
            larchname = "%s.lar" % name
            for dirname in _sys.path:
                if larchname in os.listdir(dirname):
                    isLarch = True
                    modname = os.path.abspath(os.path.join(dirname,larchname))
                    # print("Found larch module:" , modname)

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
                        if self.error:
                            break
                        
                    thismod = _sys.modules[name]               
                    symtable._set_local_mod(saveGroups)
            if self.error:
                _sys.modules.pop(name)
                thismod = None
                return

            # or, if not a larch module, load as a regular python module
            if not isLarch:
                try:
                    __import__(name)
                    thismod = sys.modules[name]
                except:
                    # print(" import error ", name, sys.exc_info())
                    self.addException(None,msg='Import Error',
                                      py_exc=sys.exc_info())
                    
                    
                    return
               
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
    # end of import_module

    # not yet implemented:
    def doExceptHandler(self,node): # ('type', 'name', 'body')
        # print("except handler")
        return None
    
    def doTryExcept(self,node):    # ('body', 'handlers', 'orelse') 
        for n in node.body:
            self.interp(n)
            if self.error:
                e_type,e_value,e_tback = sys.exc_info()
                for h in node.handlers:
                    handler_type = self.interp(h.type)
                    if handler_type is None or isinstance(e_type(),handler_type):
                        self.error = []
                        self._NodeAssign(h.name,e_value)
                        for b in h.body: self.interp(b)
                        break

                    
    def doTryFinally(self,node):    return self.NotImplemented(node)        
    def doExec(self,node):          return self.NotImplemented(node)
    def doLambda(self,node):        return self.NotImplemented(node)
    def doClass(self,node):         return self.NotImplemented(node)
    def doGlobal(self,node):        return self.NotImplemented(node)
    def doGenerators(self,node):    return self.NotImplemented(node)
    def doYield(self,node):         return self.NotImplemented(node)
    def doDecorators(self,node):    return self.NotImplemented(node)
    def doGeneratorExp(self,node):  return self.NotImplemented(node)        

#         print('Incomplete GeneratorExp ')
#         print(ast.dump(node.elt),include_attributes=True)
#         for n in node.generators:
#             print(n)             

    def doRaise(self,node):    # ('type', 'inst', 'tback') 
        msg = "%s: %s" % (self.interp(node.type).__name__,
                          self.interp(node.inst))
        self.addException(node.type,msg=msg,
                          py_exc=sys.exc_info())
