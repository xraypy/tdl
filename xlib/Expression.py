#!/usr/bin/python
#
# M. Newville Univ of Chicago (2006)
#
# Basic support for mathematical expressions in tdl language.
# Contains three principle classes:
#  opcodes:
#     simple class with attributes for operation code types.
#  ExpressionParser:
#     parse text of expression into postfix expression stack for
#     later evaluation
#  Expression:
#     handles compiling (with ExpressionParser) and evaluation
#     of tdl expressions.
#
##########################################################################
import Num
import types

import pyparsing
from pyparsing import nums, alphas, quotedString, restOfLine
from pyparsing import Word, Combine, Optional, Literal, CaselessLiteral
from pyparsing import OneOrMore, ZeroOrMore, Forward, Or
from pyparsing import QuotedString, StringEnd

from Util import trimstring, list2array
from Util import ParseError, EvalError
from Symbol import SymbolTable, symTypes,Symbol

__version__ = '0.3.4'

# opcodes for expression parsing and evaluation
class opcodes:
    function =  "@FCN"
    arrayfunc = "@AFN"
    command =   "@CMD"
    variable =  "@VAR"
    array =     "@ARR"
    string =    "@STR"
    list =      "@LIS"
    dict =      "@DIC"
    uminus =    "@UN-"
    uplus =     "@UN+"
    subarray =  "@SUB"
    slice =     "@SLI"
    slice3 =    "@SL3"
    assign =    "@ASN"
    empty =     "@EMP"
    colon =     "@COL"
    symbol =    "@SYM"
    eof =       "@EOF"
    prefix =    "@"


def make_array(n,w):
    a = []
    for i in range(n): a.append(w.pop())
    a.reverse()
    return a

def get_slice(work,use_step=False):
    " generate array slice for sytnax like a[1:4], a[2:10:2], etc"
    step  = 1
    if use_step:
        step  = work.pop()
        if step == opcodes.empty:  step = 1

    x2 = work.pop()
    if x2 == opcodes.empty:  x2 = None
    if x2 is not None:       x2 = int(x2)

    x1 = work.pop()
    if x1 == opcodes.empty:  x1 = 0
    x1 = int(x1)
    return slice(x1,x2,step)

def take_subarray(val, ndim, work):
    # takes subarrays / slices of lists and numeric arrays,
    # including syntax like: x[:2], x[3:], x[2:5,3], x[2:40:2,4:8]
    elems = []

    for i in range(ndim):
        e = work.pop()
        if type(e) == types.ComplexType: e = int(e.real)
        if type(e) == types.FloatType:   e = int(e)
        elems.append(e)
    elems.reverse()
    if type(val) == Num.ArrayType:
        if len(elems) > len(val.shape): elems = elems[:len(val.shape)]
        val = val[tuple(elems)]
    else:
        for e in elems: val = val[e]

    return val

########################################################
##
## Parser
##
_point  = Literal(".")
_equal  = Literal("=")
_comma  = Literal(",")
_colon  = Literal(":")
_percent= Literal("%")

_lpar  , _rpar   = (Literal("("), Literal(")"))
_lbrace, _rbrace = (Literal("{"), Literal("}"))
_lbrack, _rbrack = (Literal("["), Literal("]"))

_op_add = Literal("+")  | Literal("-")
_op_una = Literal('+')  | Literal('-') |  \
          Literal('!')  | CaselessLiteral('not ')
_op_mul = Literal("*")  | Literal("/")  | Literal("%")
_op_exp = Literal("^")  | Literal("**")
_op_eq  = Literal("==") | Literal("!=") | \
          Literal(">=") | Literal("<=") | \
          Literal('&&') | Literal("||") | \
          Literal("<")  | Literal(">")
_op_and = CaselessLiteral("and")
_op_or  = CaselessLiteral("or")

unary_operators = {'!':'not','not':'not',
                   '-':opcodes.uminus,'+':opcodes.uplus}


class ExpressionParser:
    """
    Expression Parser: convert a mathematical expression to a
    list of postfix tokens for easy evaluation.  The generated
    token list can be evaluated using the Evaluator class,
    which uses a SymbolTable for variable / function lookup.

    See test_expr.py for parsing & evaluation tests

    Syntax notes:
        simple mathematical operations (+,-,*,/) are supported, with
        either '**' or '^' for exponentiation, unary minus, and '%' for
        modulo.

        comparison operators are  '==', '>','<','>=', '<=', 'or', and 'and'
        symbol names

    based on the demonstration program fourFn.py from pyparsing.
    """
    def __init__(self,**kw):
        self.debug = 0
        expr, term        = Forward(), Forward()
        _slice,arg_list   = Forward(), Forward()
        lit_list,lit_dict = Forward(), Forward()

        # vname = name part (name of simple variable , group, or method)
        # name  = full name ( group:var.member, var.member, group:var,  etc)
        alphasx = alphas+"_"+"&"+"$"+"@"
        vname = Word(alphasx, alphasx+nums)
        name  = Combine(ZeroOrMore(vname + _point) +  vname)

        fnum  = (Combine(Word(nums,nums) +
                         Optional(_point + Optional(Word(nums,nums))) +
                         Optional(CaselessLiteral("e") + Word("+-"+nums, nums)) +
                         Optional(CaselessLiteral("j")))
                 ) | ( Combine((_point + Word(nums,nums)) +
                               Optional(CaselessLiteral("e") + Word("+-"+nums, nums)) +
                               Optional(CaselessLiteral("j"))) )

        _str  = ( QuotedString("'''", multiline=True).setParseAction(self.pushString_SM) |
                  QuotedString('"""', multiline=True).setParseAction(self.pushString_DM) |
                  QuotedString("'").setParseAction(self.pushString_SS) |
                  QuotedString('"').setParseAction(self.pushString_DS) )

        _num  = fnum.setParseAction(self.pushNum)
        _sym  = (name + Optional(_lpar + arg_list + _rpar)
                 + ZeroOrMore(_lbrack + _slice + _rbrack)).setParseAction(self.pushSymbol)
        _expr = ((_lpar  + expr + _rpar) +
                 ZeroOrMore(_lbrack + _slice + _rbrack).setParseAction(self.pushSubArray))
        _list = (_lbrack + lit_list + _rbrack).setParseAction(self.pushList)
        _dict = (_lbrace + lit_dict + _rbrace).setParseAction(self.pushDict)

        _emptylist = (_lbrack + _rbrack).setParseAction(self.pushEmptyList)
        _emptydict = (_lbrace + _rbrace).setParseAction(self.pushEmptyDict)

        atom  = _str | _num | _sym  | _expr | _list | _dict | _emptylist | _emptydict

        # define exponentiation as "atom [ ^ expr ]..." instead of "atom [ ^ atom ]...",
        # to get right associative
        atom = atom + ZeroOrMore( (_op_exp + term).setParseAction(self.pushBinOp) )
        term << atom
        # put Unary operations next:
        term = (OneOrMore(_op_una) + term).setParseAction(self.pushUnary) | atom

        # add other operators in order of precedence.
        for op in (_op_mul,_op_add,_op_eq,_op_and,_op_or):
            term = term + ZeroOrMore( (op + term).setParseAction(self.pushBinOp) )

        # argument list member, either simple expression or keyword/val assignment
        arg  = (vname + _equal + expr).setParseAction(self.pushKeyValArg) | \
               (expr.setParseAction(self.pushSimpleArg))

        dict_arg = (quotedString + _colon.suppress() + expr).setParseAction(self.pushDictArg)

        expr    << term
        lit_list<< (expr + ZeroOrMore(_comma.suppress() + expr)      ).setParseAction(self.CountArgs)
        arg_list<< (Optional(arg) + ZeroOrMore(_comma.suppress()+arg)).setParseAction(self.CountArgs)
        lit_dict<< (dict_arg + ZeroOrMore(_comma.suppress()+dict_arg))

        # slicing still seems a little hackish ....
        short_slice = (Optional(expr).setParseAction(self.pushSlice1) + \
                       Optional(_colon).setParseAction(self.pushColon) + \
                       Optional(expr).setParseAction(self.pushSlice2) )

        long_slice  = short_slice + Optional((_colon + Optional(expr)).setParseAction(self.pushSlice3))
        slice_item  = (short_slice ^ long_slice)
        _slice << (slice_item + ZeroOrMore((_comma + slice_item)))

        self.expr = expr + restOfLine.setParseAction(self.pushRestOfLine) + StringEnd()
        self.expr.streamline()

    def compile(self,s,reverse=False):
        self.exprStack = []
        self.argcount  = 0
        self.dictcount = 0
        s = s.strip()
        try:
            self.expr.parseString(s)
        except (ParseError,pyparsing.ParseException):
            raise ParseError, s
        # reversing the stack is the default...
        # set reverse=True to NOT reverse the stack
        if not reverse: self.exprStack.reverse()
        if self.debug>=8:
            print ' ExprParse: compile: ', s, ' -> ', self.exprStack
        return self.exprStack

    def CountArgs(self, s, loc, toks ):
        if self.debug>=16:    print ' ExprParse: CountArgs ', len(toks)
        self.argcount = len(toks)

    def pushSubArray(self, s, loc, toks ):
        n = 0
        if len(toks)>1:
            for i in toks[1:]:
                if i in (']',','): n = n + 1
        if n>0:  self.pushOp(opcodes.subarray, count=n)

    def pushSlice1(self, s, loc, toks ):
        if len(toks)==0: self.exprStack.append(opcodes.empty)
        return []

    def pushSlice2(self, s, loc, toks ):
        # print 'push Slice 2 ', s, toks, self.exprStack
        p1  = self.exprStack.pop()
        if len(toks)==0:
            if p1 != opcodes.colon:
                self.exprStack.append(p1)
            else:
                self.exprStack.append(opcodes.empty)
                self.exprStack.append(opcodes.slice)
        else:
            if (len(self.exprStack)>0 and self.exprStack[-1]==opcodes.colon):
                self.exprStack.pop()
            self.exprStack.append(p1)
            self.exprStack.append(opcodes.slice)

    def pushSlice3(self, s, loc, toks ):
        stride = self.exprStack.pop()
        xtok = self.exprStack.pop()
        self.exprStack.append(stride)
        if xtok != opcodes.slice:
            msg = "syntax error: %s.  " % s
            msg = "%s\n   bad slice syntax " % (msg)
            raise ParseError, msg
        self.exprStack.append(opcodes.slice3)

    def pushColon(self, s, loc, toks ):
        if len(toks)>0: self.exprStack.append(opcodes.colon)

    def pushRestOfLine(self, s, loc, toks ):
        if len(toks)>0 and len(toks[0]) > 0:
            msg = "syntax error: %s.  " % s
            msg = "%s\n   unrecognized tokens: %s" % (msg,toks[0])
            raise ParseError, msg

    def pushSymbol(self, s, loc, toks ):
        """ push symbols for variable and function names:
        this is the most challenging push* routine !!!
        """
        t = []
        t.append(toks[0])
        op = opcodes.array
        n,na  = 0,None
        if self.debug>=16:
            print ' ExprParse: pushSymbol ', toks, self.argcount, self.exprStack
        if len(toks)>1:
            paren_level = 0
            if toks[1] == '(': op = opcodes.function
            for i in toks[1:]:
                if i == '(' and paren_level==0:
                    paren_level = 1
                    # op = opcodes.function
                    na = self.argcount
                    self.argcount  = 0
                elif i == ')':
                    paren_level = 0
                elif i in (']',','):
                    n = n + 1
                elif i != '[' and paren_level==0:
                    t.append(i)
        if self.debug>=8:
            print ' ExprParse: pushSymbol ', op, na, n, t

        if op == opcodes.function and n==0:  # regular function call
            return self.pushFirst(s,loc,t,op=op,count=na)
        elif op == opcodes.function and n>0: # "array function"
            return self.pushFirst(s,loc,t,op=opcodes.arrayfunc,count=na,count2=n)

        elif op == opcodes.array and n==0:   # scalar variable / variable name
            return self.pushFirst(s,loc,t,op=opcodes.variable)
        else:                                # general variable slice case
            return self.pushFirst(s,loc,t,op=op,count=n)

    def pushKeyValArg(self, s, loc, toks ):
        self.pushFirst(s,loc,toks,op=opcodes.assign)
        return toks[0]

    def pushDictArg(self, s, loc, toks ):
        t = [trimstring(toks[0]), toks[1]]
        self.pushFirst(s,loc,t,op=opcodes.assign)
        self.dictcount = self.dictcount + 1

    def pushSimpleArg(self, s, loc, toks ):
        return toks[0]

    def pushUnary(self, s, loc, toks ):
        return self.exprStack.append(unary_operators[toks[0].strip().lower()])

    def pushNum(self, s, loc, toks ):
        """ push literal number, coerce to complex or float """
        cast = float
        if toks[0].endswith('j'): cast = complex
        return self.pushFirst(s,loc,[cast(toks[0])])

    def pushString(self, s, loc, toks ,format='"%s"'):
        t= [format % i for i in toks]  # t= [format % trimstring(i) for i in toks]
        return self.pushFirst(s,loc,t,op=opcodes.string)

    def pushString_SM(self, s, loc, toks ):
        self.pushString(s,loc,toks, format="'''%s'''")

    def pushString_DM(self, s, loc, toks ):
        self.pushString(s,loc,toks, format='"""%s"""')

    def pushString_SS(self, s, loc, toks ):
        self.pushString(s,loc,toks, format="'%s'")

    def pushString_DS(self, s, loc, toks ):
        self.pushString(s,loc,toks, format='"%s"')

    def pushList(self, s, loc, toks ):
        self.pushOp(opcodes.list, self.argcount,reset=True)

    def pushDict(self, s, loc, toks ):
        self.pushOp(opcodes.dict, self.dictcount)
        self.dictcount = 0

    def pushEmptyList(self, s, loc, toks):
        self.pushOp(opcodes.list, 0, reset=True)

    def pushEmptyDict(self, s, loc, toks):
        self.pushOp(opcodes.list, 0)
        self.dictcount = 0

    def pushOp(self,op,count,count2=None,reset=False):
        if count  is not None: self.exprStack.append(count)
        if count2 is not None: self.exprStack.append(count2)
        if op:     self.exprStack.append(op)
        if reset:  self.argcount=0

    def pushBinOp(self, s, loc, toks):
        return self.pushFirst(s,loc, toks)
    
    def pushFirst(self, s, loc, toks,op=None,count=None,count2=None,reset=False,**kw):
        if self.debug>=32:
            print ' ExprParse: pushFirst ', toks, op, count, count2, reset, self.exprStack
        if toks:
            self.exprStack.append(toks[0])
        if op:    self.pushOp(op,count,count2=count2,reset=reset)
        return toks

    def pushLast(self, s, loc, toks,op=None,count=None,count2=None,reset=False):
        if self.debug>=32:
            print ' ExprParse: pushLast ', toks, op, count, count2, reset, self.exprStack
        if toks:  self.exprStack.insert(0,toks[0])
        if op:    self.pushOp(op,count,count2=count2,reset=reset)
        return toks


####################################################
class Expression:
    """
    Compile and Evaluate tdl expressions (typically right-hand-side of an equation)

    The basic interface uses
    >>> e = Expression()

    where symbolTable contains the SymbolTable instance for function/variable names

    >>> stack = e.compile('2 / 5')
    >>> print stack
    ['/', 5.0, 2.0]
    >>> print e.eval(stack)
    0.4

    The .compile function compiles the expression to a postfix stack representation of
    the expression. The .eval function evaluates this stack to an actual value.
    You can also give the .eval() function an expression:

    >>> print e.eval('3/5')
    0.6
    """
    def __init__(self,symbolTable=None,run_procedure=None,debug=0):
        self.symbolTable = symbolTable
        if self.symbolTable is None: self.symbolTable = SymbolTable()

        self.debug   = debug
        self.run_procedure = run_procedure
        self.Parser  = ExpressionParser()
        self.Parser.debug = debug
        self.text    = ''


    def compile(self,expr='',reverse=False):
        if len(expr)<1: return None
        self.text = expr
        r =self.Parser.compile(expr,reverse=reverse)
        if self.debug>=4:
            print 'expr compile = %s ' % expr
            print 'expr compile -> ', r

        return r

    def set_debug(self,n):
        self.debug = n
        self.Parser.debug = n

    def check_retval(self,val):
        # checks that return value is ok by raise exception for many error conditions
        if type(val) == types.NotImplementedType: self.raise_error('Not a Number')
        if type(val) == Num.ArrayType:
            if True in Num.isnan(val): self.raise_error('Not a Number')
            if True in Num.isinf(val): self.raise_error('Infinity')
        if (type(val) in Num.typeDict.values() or \
            type(val) in [types.FloatType,types.ComplexType]):
            if Num.isnan(val):  self.raise_error('Math Error (undefined number)')
            if Num.isinf(val):  self.raise_error('Infinity')
        if type(val) == Num.ArrayType and len(val) == 1: val = val[0]
        return val

    def raise_error(self,msg):
        if len(self.text)>0: msg =  "%s at line:\n  '%s'" % (msg,self.text)
        raise EvalError, msg

    def get_symbol(self,symbol,vtype='variable',expr=''):
        if expr != '': self.text = expr
        # print 'get_symbol ', symbol, type(symbol), len(symbol), vtype
        x = None
        if vtype in ('variable','defvar'):
            x = self.symbolTable.getSymbol(symbol)
        elif vtype == 'function':
            x = self.symbolTable.getFunction(symbol)

        if x is None:
            self.raise_error( '%s "%s" not found' % (vtype,symbol) )
        return x


    def eval(self,stack=None,expr=''):
        "evaluate expression stack, compiled with Parser.parse()"

        if type(stack)==types.StringType and expr=='':  stack,expr = None,stack

        if stack is None: stack = self.compile(expr)
        if stack is None or len(stack)<1 or type(stack) != types.ListType:
            self.raise_error('cannot evaluate expression %s, %s' % (stack,expr))
        if expr != '': self.text = expr
        work = []
        code = stack[:]

        if self.debug>=8:        print ' evaluate ', expr, '\n -> ', code

        # print ' evaluate ', expr, '\n -> ', code
        # if opcodes.array in code or opcodes.subarray in code:    print 'Array Code: ', code

        while len(code)> 0:
            val = tok = code.pop()
            if self.debug>=64: print 'TOK ', tok
            if tok is None:
                self.raise_error( 'evaluation error (unrecognized expression)')
            # numeric types get pushed immediately
            if type(tok) == types.StringType:
                if tok.startswith(opcodes.prefix):
                    if tok == opcodes.colon:
                        continue
                    elif tok == opcodes.empty:
                        pass
                    elif tok == opcodes.symbol:   # special case for assignments...
                        return work
                    elif tok == opcodes.variable: # simple variable reference
                        
                        nam = work.pop()
                        sym = self.get_symbol(nam)
                        try:
                            if sym.type == symTypes.defvar: sym.value = self.eval(sym.code)
                            val = sym.value
                        except AttributeError:
                            val = sym
                    elif tok == opcodes.array: # for simple array slices: variable[slice]
                        ndim = work.pop()
                        val  = work.pop()        
                        sym = self.get_symbol(val)
                        try:
                            if sym.type == symTypes.defvar:
                                # re-evaluate defined variables here
                                sym.value = self.eval(sym.code)
                            val = sym.value
                        except AttributeError:
                            val = sym

                        val  = take_subarray(val,ndim,work)

                    elif tok == opcodes.subarray: # subarray used for (expr)[slice] and fcn(...)[slice]
                        ndim = work.pop()
                        tmp = [work.pop() for i in range(ndim)]
                        tmp.reverse()
                        val  = take_subarray(work.pop(),ndim,tmp)

                    elif tok in (opcodes.function,opcodes.arrayfunc,opcodes.command):
                        if tok == opcodes.arrayfunc:  ndim = work.pop()
                        nargs = work.pop()
                        fname = work.pop()
                        fcn   = self.get_symbol(fname,vtype='function')
                        arr = []; kws = {}; elems = []
                        if tok == opcodes.arrayfunc:
                            for i in range(ndim):  elems.append(work.pop())
                            elems.reverse()

                        for i in range(nargs):
                            x = work.pop()
                            if x==opcodes.assign:
                                y = work.pop()
                                kws[y] = work.pop()
                            else:
                                arr.append(x)
                        arr.reverse()
                        # print 'PYFUNC :: ',fcn, isinstance(fcn,Symbol),callable(fcn)
                        if isinstance(fcn,Symbol):
                            if fcn.type==symTypes.defpro and self.run_procedure:
                                val = self.run_procedure(fcn,args=arr,kws=kws)
                            elif fcn.type==symTypes.pyfunc:
                                # print 'PyFunc Call ', arr,  kws
                                val = fcn.call(*arr,**kws)
                            else:
                                self.raise_error('evaluation error for function %s' % fname)
                                
                        elif callable(fcn):
                            try:
                                val = fcn(*arr,**kws)
                            except:
                                self.raise_error('evaluation error for function %s' % fname)
                        if tok == opcodes.arrayfunc:
                            for i in elems:
                                if type(i) == types.FloatType: i = int(i)
                                val = val[i]
                        if tok == opcodes.command:
                            val = fcn.cmdout(val,**kws)

                    elif tok in (opcodes.uminus,opcodes.uplus):
                        val = work.pop()
                        if tok == opcodes.uminus: val = - val

                    elif tok == opcodes.string:
                        val = trimstring(work.pop())

                    elif tok == opcodes.list:
                        nx  = work.pop()
                        val = make_array(nx,work)

                    elif tok == opcodes.dict:
                        nx = work.pop()
                        kws = {}
                        for i in range(nx):
                            x = work.pop()
                            if x !=opcodes.assign:
                                self.raise_error('unexpected dictionary item : %s' % repr(x))
                            y = work.pop()
                            kws[y] = work.pop()
                        val = kws

                    elif tok in (opcodes.slice,opcodes.slice3):
                        val = get_slice(work,use_step=(tok == opcodes.slice3))
                    elif tok != opcodes.assign:
                        self.raise_error('unknown token : %s in "%s"' % (tok , expr))

                elif tok == '%':
                    # may mean modulo or string-format...
                    x = work.pop() ;  y = work.pop()
                    if type(y) == types.StringType  and \
                       type(x) in (types.ListType, Num.ArrayType): x = tuple(x)
                    val = y %  x

                elif tok in ('!', 'not'): x = work.pop() ; val = not x
                elif tok in ('^','**'):   x = work.pop() ; val = work.pop() ** x
                elif tok == '!':   x = work.pop() ; val = not x
                elif tok == '+':   x = work.pop() ; val = work.pop() +  x
                elif tok == '-':   x = work.pop() ; val = work.pop() -  x
                elif tok == '*':   x = work.pop() ; val = work.pop() *  x
                elif tok == '/':   x = work.pop() ; val = work.pop() /  x
                elif tok == '==':  x = work.pop() ; val = work.pop() == x
                elif tok == '!=':  x = work.pop() ; val = work.pop() != x
                elif tok == '>':   x = work.pop() ; val = work.pop() >  x
                elif tok == '<':   x = work.pop() ; val = work.pop() <  x
                elif tok == '>=':  x = work.pop() ; val = work.pop() >= x
                elif tok == '<=':  x = work.pop() ; val = work.pop() <= x
                elif tok == 'or':  x = work.pop() ; val = work.pop() or x
                elif tok == 'and': x = work.pop() ; val = work.pop() and x
                #
            if type(val)==types.TupleType: val = list(val)
            work.append(val)
            if self.debug >=32: print ' work ', work
        #
        # print 'WORK ', work, type(work), len(work)

        if type(work) == types.ListType:
            work = list2array([self.check_retval(i) for i in work])
        if len(work)==1:
            work = work[0]
        if type(work) == types.StringType:
            work = trimstring(work)
        return work

if __name__ == '__main__':
    p = Expression(debug=1)
    s = p.symbolTable
    s.setVariable('a',Num.arange(30.))
    s.addGroup('g1')
    s.setVariable('g1.a',3.2)
    s.setVariable('g1.b',99.0)
    s.setVariable('b',Num.arange(15.))
    s.setVariable('s1','A long StrinG')
    s.setFunction('sqrt',Num.sqrt)
    s.setVariable('format1',' %s = %f ')
    s.setVariable('dlist',['abcdefghijklmnop',12,'xyz'])
    s.setVariable('adict',{'key1': 'abcdefghijklmnop','key2':3.1})
    s.setVariable('x',2.2)
    s.setVariable('ix',4)
    s.setVariable('st','a long string here')

    t = ('a', 'a[2]', 'a[:7]', 'a[4:10]', 'b[9:]')
    t = ('a', 'a[2]',
         '""" three q"""',
         "''' three q  '''",
         "'strform %s ' % 's'",
         "'strform %s ' % ('s')",
         "' %s = %f ' % ['a',3.3]",
         " format1  % ['a',3.3]")
    #" format1  % dlist",         )
    st = (' (x+1)', 'sqrt((x+1)) ' ,
         'sqrt((a+1)/4)[3]' , 'sqrt(x)' ,
         'st[2:8]',
         'a[2]',
         'adict["key1"][1:4]',
         'dlist[0,1:4]',
         '[3, 4, 5]',
         )
    t = ('b', 'g1.a','sqrt(g1.b)','sqrt(b)','s1', 's1.upper()','(x+1)')
    t = ('b', 'x', '(x+1)','2*3/4')
    for i in t:
        print '========================\n< ', i , ' > '
        x = p.compile(i)
        print '--> ', x
        y = p.eval(x)
        print '=> ', y
