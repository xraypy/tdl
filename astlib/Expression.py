import os, sys
from compiler import ast, parse
import operator

import Num
import Symbol

BinOps = {"is": operator.is_, "is not": operator.is_not, "in": operator.contains, 
          "==": operator.eq , "!=": operator.ne,
          "<" : operator.lt , ">" : operator.gt, "<=": operator.le , ">=": operator.ge}

AugAssignOps = {'+=': operator.add, '-=': operator.sub,
                '*=': operator.mul, '/=': operator.div,
                '%=': operator.mod, '**=': operator.pow,
                '>>=': operator.rshift, '<<=': operator.lshift,
                '&=':  operator.and_, '|=':  operator.or_ }
                  

class opcodes:
    delim1,delim2   = "@;",";@"
    _continue ="%scontinue%s"%(delim1,delim2)
    _pass     ="%spass%s"   %(delim1,delim2)
    _break    ="%sbreak%s"  %(delim1,delim2)
    _print    ="%sprint%s"  %(delim1,delim2)
    _return   ="%sreturn%s" %(delim1,delim2)
    _raise    ="%sraise%s"  %(delim1,delim2)
    _assert   ="%sassert%s" %(delim1,delim2)
    _eof      ="%seof%s"    %(delim1,delim2)
    minlen    = len(delim1+delim2)

def valid_opcode(t):
    return (t.startswith(opcodes.delim1) and
            t.endswith(opcodes.delim2) and
            len(t)>opcodes.minlen)
   
def isGroup(obj):  return isinstance(obj,Symbol.Group)
def isSymbol(obj): return isinstance(obj,Symbol.Symbol)

class EvalError(Exception):
    def __init__(self,error,descr = None,node = None):
        self.error = error
        self.descr = descr
    def __repr__(self):
        return "%s: %s" % (self.error, self.descr)
    __str__ = __repr__


class Compiler:
    """ Statement compiler and interpreter.
    This module compiles expressions and statements to AST representation,
    relying on python's syntax and the compiler module, and then executes these
    representations using a SymbolTable for named objects (variable, functions).

    This is a restricted version of Python, and does not include any block structures,
    or flow control (loops, conditionals), etc.  Those bits are reserved for a calling
    class (such as Interpretor).

    In addition, this module supports a sub-set of python syntax. The following syntax is not supported:
        Ellipses [...] for array / list access
        List comprehension, Generators, Yield, Decorators,
        Class, Import, From, Exec
        Global, Function, Lambda
        If, While, TryExcept, TryFinally, ...        
    
    """
    def __init__(self):
        self.symbolTable = Symbol.SymbolTable()
        self.__methods = {
            ast.Module: self.doModule,
            ast.Stmt: self.doStmt,
            ast.Expression: self.doExpression,
            ast.Discard: self.doDiscard,
            ast.Add: self.doAdd,
            ast.Sub: self.doSub,
            ast.Mul: self.doMul,
            ast.Div: self.doDiv,
            ast.FloorDiv: self.doFloorDiv,
            ast.Power: self.doPower,
            ast.Mod: self.doMod,
            ast.LeftShift: self.doLeftShift,
            ast.RightShift: self.doRightShift,
            ast.And: self.doAnd,
            ast.Or: self.doOr,
            ast.Bitand: self.doBitand,
            ast.Bitor: self.doBitor,
            ast.Bitxor: self.doBitxor,
            ast.Not: self.doNot,
            ast.Invert: self.doInvert,
            ast.UnaryAdd: self.doUnaryAdd,
            ast.UnarySub: self.doUnarySub,
            ast.Backquote: self.doBackquote,
            ast.Compare: self.doCompare,
            ast.Const: self.doConst,
            ast.List: self.doList,
            ast.Tuple: self.doTuple,
            ast.Dict: self.doDict,
            ast.Keyword: self.doKeyword,
            ast.Getattr: self.doGetattr,
            ast.Name: self.doName,
            ast.Sliceobj: self.doSliceobj,
            ast.Slice: self.doSlice,
            ast.Subscript: self.doSubscript,
            ast.CallFunc: self.doCallFunc,
            ast.AssAttr: self.doAssAttr,
            ast.AssList: self.doAssList,
            ast.AssName: self.doAssName,
            ast.AssTuple: self.doAssTuple,
            ast.Assign: self.doAssign,
            ast.AugAssign: self.doAugAssign,
            ast.Continue: self.doContinue,
            ast.Pass: self.doPass,
            ast.Return: self.doReturn,
            ast.Break: self.doBreak,
            ast.Raise: self.doRaise,
            ast.Print: self.doPrint,
            ast.Printnl: self.doPrintnl,

            }
        
    def NotImplemented(self,node,**kw):
        raise EvalError, "syntax error:  '%s' not supported" % (node.__class__.__name__)
    
    # main entry point for Ast node evaluation
    #  compile:  string statement -> ast
    #  interp :  ast -> result
    #  eval   :  string statement -> result = interp(compile(statement))
    def compile(self,expr,mode='single'):
        """compile statement/expression to Ast representation
        mode = single  to compile a single interactive command
             = eval    to compile an expression
             = exec    to compile a module
        """
        if mode not in ('single','eval','expr'): mode='single'
        return parse(expr,mode)
    
    def interp(self, node,**kw):
        """executes compiled Ast representation for an expression"""
        if node is None: return None
        # print 'Expr Interp ', node
        return self.__methods.get(node.__class__,self.NotImplemented)(node,**kw)

    def eval(self,expr,mode='single'):
        """evaluates a single statement"""
        return self.interp(self.compile(expr,mode=mode))

    # outer nodes: Module, Statement, Expression
    def doModule(self,node,**kw):
        'Module node, contains doc and node (a statement)'
        try:
            if node.doc is not None and len(node.node.nodes)==0:
                return node.doc
        except:
            pass
        return self.doStmt(node.node,**kw)


    def doStmt(self,node,**kw):
        'Statement node: a set of child nodes'
        o = None
        for child in node.nodes: o = self.interp(child,**kw)
        return o
            
    def doExpression(self, node, **kw):
        'Expression node: '
        o = None
        for child in node.getChildNodes():  o = self.interp(child, **kw)
        return o

    def doDiscard(self, node, **kw):
        'Discard node'
        return self.interp(node.expr,**kw)
    
    # binary ops, all having a right and left node
    def doAdd(self, node, **kw):
        'Add two nodes'
        return self.interp(node.left) + self.interp(node.right)
    
    def doSub(self,node,**kw):
        'Subtract two nodes'
        return self.interp(node.left) - self.interp(node.right)
    
    def doMul(self,node,**kw):
        'Multiply two nodes'        
        return self.interp(node.left) * self.interp(node.right)
    
    def doDiv(self,node,**kw):
        'Divide two nodes'        
        return operator.truediv(self.interp(node.left), self.interp(node.right))

    def doFloorDiv(self,node,**kw):
        'FloorDivision two nodes'        
        return self.interp(node.left) // self.interp(node.right)

    def doPower(self,node,**kw):
        'Power: exponentiate two nodes'        
        return self.interp(node.left) ** self.interp(node.right)
    
    def doMod(self,node,**kw):
        'Modulo nodes: also used for string formatting!'
        return self.interp(node.left) % tuple(self.interp(node.right))

    def doLeftShift(self,node,**kw):
        'LeftShift two nodes'
        return self.interp(node.left) << self.interp(node.right)

    def doRightShift(self,node,**kw):
        'RightShift two nodes'
        return self.interp(node.left) >> self.interp(node.right)

    def doAnd(self, node, **kw):
        'And a set of nodes'
        for arg in node.nodes:
            if not self.interp(arg): return False
        return True

    def doOr(self, node, **kw):
        'Or a set of nodes'
        for arg in node.nodes:
            if self.interp(arg): return True
        return False

    def doBitand(self, node, **kw):
        'Bitand a set of nodes'
        return reduce(lambda x,y: x&y, [self.interp(n) for n in node.nodes])

    def doBitor(self, node, **kw):
        'Bitor a set of nodes'
        return reduce(lambda x,y: x|y, [self.interp(n) for n in node.nodes])

    def doBitxor(self, node, **kw):
        'Bitxor a set of nodes'
        return reduce(lambda x,y: x^y, [self.interp(n) for n in node.nodes])

    def doNot(self,node,*kw):
        'Not a node'
        return not self.interp(node.expr)

    def doInvert(self, node, **kw):
        'Invert a node'
        return ~self.interp(node.expr)

    def doUnaryAdd(self,node,*kw):
        'UnaryAdd a node'
        return +self.interp(node.expr)

    def doUnarySub(self,node,*kw):
        'UnarySub a node'
        return -self.interp(node.expr)
    
    def doBackquote(self, node, **kw):
        'Backquote node'
        return repr(self.interp(node.expr))

    def doCompare(self, node, **kw):
        'Compare an expression with a set of node (ops,expressions)'
        val = self.interp(node.expr)
        for op, subnode in node.ops:
            args = [val, self.interp(subnode)]
            if op == 'in':
                args.reverse()  # ie, 'x in y' == contains(y,x)
            if not BinOps[op](*args): return False
        return True

    def doConst(self, node, **kw):
        'Const node: return the constant value'
        return node.value

    def doList(self, node, **kw):
        'List node: create a list from subnodes'
        return [self.interp(s) for s in node.nodes]

    def doTuple(self, node, **kw):
        'Tuple node: create a tuple from subnodes'
        return tuple(self.interp(s) for s in node.nodes)

    def doDict(self, node, **kw):
        """Dict node: build a dictionary from items in node:
        enforce keywords being strings here"""
        return dict([(str(self.interp(k)),self.interp(v)) for k,v in node.items])

    def doKeyword(self,node,**kw):
        'Keyword node (for function arguments): return key, val'
        return node.name, self.interp(node.expr)

    # call a function
    def doCallFunc(self, node, **kw):
        'CallFunc node'
        fcn = self.interp(node.node)
        kws = self.interp(node.dstar_args) or {}
        args= []
        for arg in node.args:
            if isinstance(arg, ast.Keyword):
                k,v = self.interp(arg)
                if k in kws:
                    raise TypeError, \
                          "%s() got multiple values for keyword argument '%s'" % (fcn,k)
                kws[k] = v
            else:
                args.append(self.interp(arg))
        args.extend(node.star_args or [])
        return fcn(*args,**kws)

    # Name resolution
    def doGetattr(self,node,for_assign=False, **kw):
        if for_assign and not hasatter(node, 'flags'): node.flags = 'OP_ASSIGN'
        obj = self.interp(node.expr)
        if isGroup(obj):
            sym = obj.get(node.attrname)
            if isSymbol(sym): sym = sym.value
        else:
            sym = getattr(obj,node.attrname)            
        return sym

    def doName(self,node, getvalue=True,  for_assign=False, **kw):

        sym = self.symbolTable.getSymbol(node.name)
        if sym is None and for_assign: # may be that we have to create a group here
            sym = self.symbolTable.addGroup(node.name)
        if getvalue and isSymbol(sym):  sym =  sym.value
        return sym
            
    def doSliceobj(self, node, **kw):
        'Sliceobj node: create a slice from subnodes'
        return slice(*[self.interp(s) for s in node.nodes])

        
    def doSlice(self,node,**kw):
        'Slice node'
        if node.flags == "OP_APPLY":
            obj = self.interp(node.expr)
            return  obj[self.interp(node.lower):self.interp(node.upper)]
        elif node.flags == "OP_ASSIGN":
            obj   = self._get_AssignObject(node.expr)            
            lo,hi = self.interp(node.lower),self.interp(node.upper)
            return obj, lo, hi

    def doSubscript(self,node, for_assign=False, **kw):
        'Subscript node'
        if for_assign and not hasatter(node, 'flags'): node.flags = 'OP_ASSIGN'        
        if node.flags == 'OP_ASSIGN': 
            obj  = self._get_AssignObject(node.expr)
            # Need to build tuple of
            if isSymbol(obj) and type(obj.value) == Num.ArrayType:
                snodes = [self.interp(s) for s in node.subs]
                return obj, tuple(snodes)
            else:
                last_sub = self.interp(node.subs.pop())
                for s in node.subs: obj = obj[self.interp(s)]
                return obj, last_sub
        else:
            obj = self.interp(node.expr)
            if type(obj) == Num.ArrayType:
                obj = obj[tuple([self.interp(s) for s in node.subs])]
            else:
                for s in node.subs: obj = obj[self.interp(s)]
            return obj
                   
    # Assignment Nodes
    #
    #  Which Ast Modules can be on the LHS, in an assignment:
    #    Name, Getattr, Slice, Subscript
    #  these must know when they are part of a LHS, and be ready
    #  to return the symbol reference, not its value as appropriate.
    # 
    #  Lots and Lots of testing needed.
    def doAssAttr(self, node, **kw):
        'AssAttr node: ????can also delete attribute name with OP_DELETE ??? '
        obj   = self._get_AssignObject(node.expr)
        if isGroup(obj):
            if not obj.has_key(node.attrname):
                self.symbolTable.createSymbolInGroup(obj,node.attrname)
            return (obj, node.attrname)
        else:            
            if hasattr(obj,node.attrname):
                return (obj, node.attrname)
            elif isSymbol(obj) and hasattr(obj.value,node.attrname):
                return (obj.value, node.attrname)
        raise EvalError, " cannot get attribute '%s' for assignment" % node.attrname

    def doAssName(self, node, **kw):
        'AssName node: can also delete symbol name with OP_DELETE !! '
        if node.flags == 'OP_DELETE':
            self.symbolTable.delSymbol(node.name)            
            return None
        return self.symbolTable.getSymbol(node.name,insert='Symbol')

    def doAssTuple(self, node, **kw):
        'AssTuple node'
        return tuple(self.doAssList(node, **kw))

    def doAssList(self, node, **kw):
        'AssList node'
        return [self.interp(s) for s in node.nodes]


    def _AssObjVal(self,obj,value,errstr=None):
        if isSymbol(obj):
            obj.value = value
        elif isGroup(obj) and isgroup(value):
            obj = value
        else:
            if errstr is None:
                errstr = " cannot assign '%s' value '%s'" % (repr(obj),repr(val))
            raise EvalError, errstr
        
    def doAssign(self, node, **kw):
        'Assign node'
        value = self.interp(node.expr)
        for s in node.nodes:
            cname = s.__class__.__name__
            # print 'Assign: ', cname
            if cname in ('AssTuple','AssList'):
                target = self.interp(i)
                if isinstance(target,(list,tuple)):
                    if len(target) != len(value):
                        raise EvalError, " need %i values to unpack" % len(target)
                    for t,v in zip(target,value): t.value = v
                else:
                    raise EvalError, " cannot assign to %s " % repr(target)

            elif cname == 'AssName':
                obj = self.doAssName(s)
                self._AssObjVal(obj,value)

            elif cname == 'AssAttr':
                obj,attr = self.doAssAttr(s)
                # print 'Assign/AssAttr:  returned obj, attr, value ' , obj, attr, value
                if isGroup(obj):
                    self.symbolTable.createSymbolInGroup(obj,attr,value=value)
                else:
                    setattr(obj,attr,value)

            elif cname == 'Subscript':
                obj,subs  = self.doSubscript(s)
                obj[subs] = value

            elif cname == 'Slice':
                obj, lo, hi = self.doSlice(s)

            else:        
                raise EvalError, " cannot assign with %s " % repr(s)
                
        
    def _get_AssignObject(self,node):
        """ return an object for assignment, as on LHS of statements"""
        node_class = node.__class__

        fcns = {ast.Name:    self.doName,
                ast.Getattr: self.doGetattr,
                ast.Subscript: self.doSubscript}

        # print '_get_AssignObject 1 ', node_class, ' :: ', node_class in fcns.keys()
        if node_class not in fcns.keys():
            raise EvalError, "unidentified object for Assignment %s " % node
        try:
            return self.__methods[node_class](node, getvalue=False,for_assign=True)
        except:
            raise EvalError, "unidentified object for Assignment %s " % node

    def doAugAssign(self, node, **kw):
        'AugAssign node: a fairly simple case??? '

        anode = node.node
        node_class = anode.__class__

        op = AugAssignOps.get(node.op,None)
        if op is None:
            raise EvalError, "unsupported Augmented Assignmet '%s'" % (node.op)

        change = self.interp(node.expr)
        curval = self.interp(anode)
        value  = op(curval,change)

        if node_class == ast.Name:
            obj =self.symbolTable.getSymbol(anode.name,insert='Symbol')
            return self._AssObjVal(obj, value)
        elif node_class in (ast.Subscript, ast.Slice):
            anode.flags = 'OP_ASSIGN'
            obj,subs  = self.doSubscript(anode)
            obj[subs] = value
               
        elif node_class == node.Getattr:
            obj = self.doAssAttr(anode)
            self._AssObjVal(obj,value)
        return None

    # these all send special signals back to the caller
    # as they indicate interruptions in execution flow
    # or (print) need special output handling

    def doContinue(self, node, **kw):   return opcodes._continue
    def doPass(self, node, **kw):       return opcodes._pass
    def doBreak(self, node, **kw):      return opcodes._break
    
    def doAssert(self, node, **kw):
        return (opcodes._assert, self.interp(node.test),self.interp(node.fail))

    def doReturn(self, node, **kw):
        return (opcodes._return, self.interp(node.value))

    def doRaise(self, node, **kw):
        return (opcodes._raise, self.interp(node.expr1),
                self.interp(node.expr2), self.interp(node.expr3) )
    
    def doPrint(self, node, **kw):
        out = [repr(self.interp(s)) for s in node.nodes]
        return (opcodes._print, self.interp(node.dest), False, out)
    
    def doPrintnl(self, node, **kw):
        out = [repr(self.interp(s)) for s in node.nodes]
        return (opcodes._print, self.interp(node.dest), True, out)


if __name__ == '__main__':
    e = Compiler()
    s = e.symbolTable

    s.setVariable('a',Num.arange(30.))
    array2d = Num.arange(36) *0.5
    array2d.shape = (9,4)
    s.addGroup('g1')
    s.setVariable('g1.a',3.2)
    s.setVariable('g1.b',99.0)
    s.setVariable('b',array2d)
    s.setVariable('s1','A long StrinG')
    s.setFunction('sqrt',Num.sqrt)
    s.setVariable('format1',' "%s" = %f ')
    s.setVariable('dlist',['abcdefghijklmnop',12,'xyz'])
    s.setVariable('adict',{'key1': 'abcdefghijklmnop','key2':3.1})
    s.setVariable('x',2.2)
    s.setVariable('ix',4)
    s.setVariable('st','a long string here')

    t = ['1.1', '0.2','0.3j','.3e3',
         '2*3/4',         '[3, 4, 5]',
         'a','g1.a',
         'a[2]', 'a[:7]', 'a[4:10]', 'b[9:]',
         '"a string"',
         "''' three q  '''",
         '""" three q"""',
         "''' three q  '''",
         "'strform %s ' % 's'",
         "'strform %s ' % ('s')",
         ' (x+1)', 'sqrt((x+1)) ' ,
         'sqrt((a+1)/4)[3]' , 
         'st[2:8]',
         'adict["key1"][1:4]',
         'dlist[0,1:4]',
         'sqrt(g1.b)','sqrt(b)','s1', 's1.upper()',
         'dlist', 'dlist[0].upper()[4:7]', 'b[3:9]/2.7',
         "'%s = %f ' % ['a',3.3]",
         "format1  % ['a',5.412]",
         "'xyz' in dlist",
         ]

    for i in t:
        i = i.strip()
        print '======\n"', i , '"  -> ', 
        x = e.compile(i)
        print x
        y = e.interp(x)
        print '=> ', y

