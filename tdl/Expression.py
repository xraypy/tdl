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
from Num import Num
import types

from pyparsing import nums, alphas, quotedString, restOfLine
from pyparsing import Word, Combine, Optional, Literal, CaselessLiteral
from pyparsing import OneOrMore, ZeroOrMore, Forward, Or
from pyparsing import QuotedString, StringEnd, ParseException

from Util import trimstring, list2array, int2bin
from Util import ParseException, EvalException
from Symbol import SymbolTable

__version__ = '0.3.2'

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

def take_subarray(val,elems):
    # takes subarrays / slices of lists and numeric arrays,
    # including syntax like: x[:2], x[3:], x[2:5,3], x[2:40:2,4:8]
    el0  = elems.pop(0)
    if type(el0) == types.FloatType: el0 = int(el0)
    if el0 != opcodes.empty: val = val[el0]

    if type(val) == Num.arraytype:
        j = 0
        for el in elems:
            j = j + 1
            if type(el) == types.SliceType:
                val = val.swapaxes(0,j)[el].swapaxes(0,j)
            elif el != opcodes.empty:
                val = val.swapaxes(0,j)[int(el)]
    else:
        for el in elems: val = val[int(el)]
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
        vname = Word(alphas+"_"+"&"+"$", alphas+nums+"_")
        name  = Combine(ZeroOrMore(vname + _point) +  vname)
        
        fnum  = (Combine(Word(nums,nums) +
                         Optional(_point + Optional(Word(nums,nums))) +
                         Optional(CaselessLiteral("e") + Word("+-"+nums, nums)) +
                         Optional(CaselessLiteral("j")))
                 ) | ( Combine((_point + Word(nums,nums)) +
                               Optional(CaselessLiteral("e") + Word("+-"+nums, nums)) +
                               Optional(CaselessLiteral("j"))) )
        
        _str  = ((QuotedString("'''", multiline=True)|
                  QuotedString('"""', multiline=True)|
                  quotedString)).setParseAction(self.pushString)

        _num  = fnum.setParseAction(self.pushNum) 
        _sym  = (name + Optional(_lpar + arg_list + _rpar)
                 + ZeroOrMore(_lbrack + _slice + _rbrack)).setParseAction(self.pushSymbol)
        _expr = ((_lpar  + expr + _rpar) +
                ZeroOrMore(_lbrack + _slice + _rbrack).setParseAction(self.pushSubArray))
        _list = (_lbrack + lit_list + _rbrack).setParseAction(self.pushList)
        _dict = (_lbrace + lit_dict + _rbrace).setParseAction(self.pushDict)

        atom  = _str | _num | _sym  | _expr | _list | _dict  

        # define exponentiation as "atom [ ^ expr ]..." instead of "atom [ ^ atom ]...",
        # to get right associative
        atom = atom + ZeroOrMore((_op_exp + term).setParseAction(self.pushFirst))
        term << atom 
        # put Unary operations next:
        term = (OneOrMore(_op_una) + term).setParseAction(self.pushUnary) | atom
        
        # add other operators in order of precedence.
        for op in (_op_mul,_op_add,_op_eq,_op_and,_op_or):
            term = term + ZeroOrMore((op + term).setParseAction(self.pushFirst))

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
        except ParseException:
            self.exprStack = []
            raise ParseException, s
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
                if i == ']':  n = n + 1
        if n>0: self.pushOp(opcodes.subarray, count=n)

    def pushSlice1(self, s, loc, toks ):
        if len(toks)==0: self.exprStack.append(opcodes.empty)
        return []

    def pushSlice2(self, s, loc, toks ):
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
            raise ParseException, msg
        self.exprStack.append(opcodes.slice3)

    def pushColon(self, s, loc, toks ):
        if len(toks)>0: self.exprStack.append(opcodes.colon)

    def pushRestOfLine(self, s, loc, toks ):
        if len(toks)>0 and len(toks[0]) > 0:
            msg = "syntax error: %s.  " % s
            msg = "%s\n   unrecognized tokens: %s" % (msg,toks[0])
            raise ParseException, msg

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
            for i in toks[1:]:
                if i == '(' and paren_level==0:
                    paren_level = 1
                    op = opcodes.function
                    na = self.argcount
                    self.argcount  = 0
                elif i == ')':
                    paren_level = 0
                elif i in (']',','):
                    n = n + 1
                elif i != '[' and paren_level==0:
                    t.append(i)
        if self.debug>=16:
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

    def pushString(self, s, loc, toks ):
        if loc == 0 and (s.startswith('"""') or s.startswith("'''")):
            t = ["'''%s'''" % trimstring(i) for i in toks]
        else:
            t = ["'%s'" % trimstring(i) for i in toks]
        return self.pushFirst(s,loc,t,op=opcodes.string)
    
    def pushList(self, s, loc, toks ):
        self.pushOp(opcodes.list, self.argcount,reset=True)

    def pushDict(self, s, loc, toks ):
        self.pushOp(opcodes.dict, self.dictcount)
        self.dictcount = 0

    def pushOp(self,op,count,count2=None,reset=False):
        if count  != None: self.exprStack.append(count)
        if count2 != None: self.exprStack.append(count2)
        if op:     self.exprStack.append(op)            
        if reset:  self.argcount=0
        
    def pushFirst(self, s, loc, toks,op=None,count=None,count2=None,reset=False):
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
        if self.symbolTable == None: self.symbolTable = SymbolTable()

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
        if type(val) in (types.FloatType,types.ComplexType):
            if Num.isnan(val):  self.raise_error('Not a Number')
            if Num.isinf(val):  self.raise_error('Infinity')
        return val
            
    def raise_error(self,msg):
        if len(self.text)>0: msg =  "%s at line:\n  '%s'" % (msg,self.text)
        raise EvalException, msg
        
    def get_symbol(self,symbol,type='variable',expr=''):
        if expr != '': self.text = expr
        if type in ('variable','defvar'):
            x = self.symbolTable.getVariable(symbol.strip())
        elif type == 'function':
            x = self.symbolTable.getFunc(symbol.strip())
        if x is None:
            self.raise_error( '%s "%s" not found' % (type,symbol) )
        try:
            return x
        except:
            self.raise_error( 'invalid %s: %s' % (type,symbol) )

    def eval(self,stack=None,expr=''):
        "evaluate expression stack, compiled with Parser.parse()"

        if type(stack)==types.StringType and expr=='':  stack,expr = None,stack

        if stack==None: stack = self.compile(expr)
        if stack==None or len(stack)<1 or type(stack) != types.ListType:
            self.raise_error('cannot evaluate expression %s, %s' % (stack,expr))
        
        if expr != '': self.text = expr
        
        work = []        
        code = stack[:]
        if self.debug>=8:        print ' evaluate ', expr, '\n -> ', code

        while len(code)> 0:
            val = tok = code.pop()
            if self.debug>=64: print 'TOK ', tok
            
            if tok==None:
                self.raise_error( 'evaluation error (unrecognized expression)')

            # numeric types get pushed immediately
            if type(tok) == types.StringType:
                if tok.startswith(opcodes.prefix):
                    if tok == opcodes.empty:
                        pass
                    ## special case for assignments...
                    elif tok == opcodes.symbol:
                        return work

                    elif tok == opcodes.variable: # simple variable reference
                        nam = work.pop()
                        sym = self.get_symbol(nam)
                        if sym.type == 'defvar': sym.value = self.eval(sym.code)
                        val = sym.value

                    elif tok == opcodes.array: # array slice / subarray
                        # handles name[slice] syntax
                        ndim  = work.pop()
                        nam   = work.pop()
                        elems = make_array(ndim,work)
                        
                        sym   = self.get_symbol(nam)
                        if sym.type == 'defvar': sym.value = self.eval(sym.code)
                        val   = take_subarray(sym.value,elems)

                    elif tok == opcodes.subarray:
                        # handles (expr)[slice] and fcn(...)[slice]
                        # syntax... that is, to take a subarray from
                        # the values in the current work array
                        ndim  = work.pop()
                        elems = make_array(ndim,work)
                        val   = take_subarray(work.pop(),elems)

                    elif tok in (opcodes.function,opcodes.arrayfunc,opcodes.command):
                        if tok == opcodes.arrayfunc:  ndim = work.pop()
                        nargs = work.pop()
                        fname = work.pop()
                        fcn   = self.get_symbol(fname,type='function')
                        
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
                        if fcn.type=='defpro' and self.run_procedure:
                            val = self.run_procedure(fcn,args=arr,kws=kws)
                        elif fcn.type=='pyfunc':
                            val = fcn(*arr,**kws)
                        else:
                            self.raise_error('evaluation error for function %s' % fname)
                        if tok == opcodes.arrayfunc:
                            for i in elems:
                                if type(i) == types.FloatType: i = int(i)
                                val = val[i]
                        if tok == opcodes.command:
                            val = fcn.__cmdout__(val)
                            
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
            work.append(val)
            if self.debug >=32: print ' work ', work
        #
        if len(work)==1:                   work = self.check_retval(work[0])
        if type(work) == types.StringType: work = trimstring(work)
        if type(work) == types.ListType:
            work = list2array([self.check_retval(i) for i in work])
        return work

if __name__ == '__main__':
    p = Expression(debug=4)
    s = p.symbolTable
    s.addVariable('a',Num.arange(30.))
    s.addVariable('b',Num.arange(15.))
    s.addVariable('format1',' %s = %f ')
    s.addVariable('dlist',['b',12])
    s.addVariable('x',2.2)
        
    t = ('a', 'a[2]', 'a[:7]', 'a[4:10]', 'b[9:]')
    t = ('a', 'a[2]',
         '""" three q"""',
         "''' three q  '''",
         "'strform %s ' % 's'",
         "'strform %s ' % ('s')",
         "' %s = %f ' % ['a',3.3]",
         " format1  % ['a',3.3]",
         " format1  % dlist",
         )
    t = (' sqrt(x+1)', 'sqrt((x+1)) ' ,
         'sqrt((a+1)/3)',)
    t = (     'sqrt((a+1)/4)[3]' , )
         
    for i in t:
        print '========================\n< ', i , ' > '
        x = p.compile(i)
        y = p.eval(x)
        print ' = ', y
