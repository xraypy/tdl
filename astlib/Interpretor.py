#!/usr/bin/env python
# M. Newville Univ of Chicago (2006)
# Evaluator:  main TDL language evaluation class
#
# --------------
# Modifications
# --------------
#
#
###########################################################################
import sys
import types
import copy

# import default library of cmds and funcs
# import Help

from Expression import Compiler, opcodes

from Symbol import SymbolTable, symTypes, Group, isGroup, isSymbol

#  import version
from Num import num_version, ArrayType
from Util import split_delim, find_unquoted_char, parens_matched
from Util import split_list, trimstring
from Util import EvalError, Command2Expr
import exceptions

class Fifo:
    def __init__(self):
        self.front,self.back = ([],[])

    def __len__(self):   return len(self.front)+len(self.back)
            
    def append(self, value): self.back.append(value)

    def pop(self):
        if not self.front:
            self.front, self.back = self.back, self.front
            self.front.reverse()
        return self.front.pop()


class Interpretor:
    """ main evaluation class for tdl language.
    parses  / compiles tdl statements and evaluates them:
    """
    EOF      = opcodes._eof
    reserved_words = ('if','else','elif', 'endif',
                      'for', 'endfor','while','endwhile',
                      'try', 'except', 'endtry', 'def', 'enddef')
                          
    __interrupts = ['pass','continue','break','return']
    def __init__(self, libs=None, interactive=True, debug=False, **kw):
        self.interactive = interactive
        self.input       = sys.stdin
        self.output      = sys.stdout
        self.debug       = debug
        # self.help        = Help.Help(tdl=self,output=self.output)

        self.Expression  = Compiler()
        self.symbolTable = self.Expression.symbolTable = SymbolTable()
        self.expr_eval   = self.Expression.eval
        self.expr_compile= self.Expression.compile
        self.expr_interp = self.Expression.interp

        self.prompt   = '... '
        self.stack = []
        self.text  = Fifo()
        self.line_buff = ''
        self.nline = 0
        self.triple_squote = False
        self.triple_dquote = False
        self.line_join     = ' '
        self.interrupt  = 0
        self.retval = None
        self.infile = '<stdin>'

        self.debug       = debug
        # self.set_debug(debug)

    def setVariable(self,var,val,**kw):
        return self.symbolTable.setSymbol(var,val,**kw)

    def getVariable(self,var):
        "return reference to a named variable"
        return self.symbolTable.getSymbol(var)

    def getVariableValue(self,var):
        "return value of a named variable"
        try:
            return self.symbolTable.getSymbol(var).value
        except AttributeError:
            return None
        return None

#     def raise_error(self,msg):
#         if len(self.text)>0: msg =  "%s at line %i:\n  '%s'" % (msg,self.nline,self.text)
#         raise EvalError, msg

    def load_file(self,fname):
        try:
            f = open(fname)
            l = f.readlines()
            f.close()
        except IOError:
            raise EvalError, 'cannot load file %s '% fname
        self.infile = fname
        self.text.append((self.EOF,-1,fname))
        self.load_statements(l,file=fname)

    def clear(self):
        " clear all statements "
        self.stack = []
        self.text  = Fifo()
        self.nline = 0

    def load_statements(self,t,file='stdin'):
        " load a list of text lines to be parsed & executed"
        if not isinstance(t,list): t = [t]
        # print 'load statement ', t
        n = 0
        self.infile = file
        while t:
            n = n+1
            s = t.pop().strip()
            self.text.append((s,n,self.infile))
            
    def execute(self,t,file='stdin'):
        " evaluate tdl statement or list of tdl statements"
        self.load_statements(t,file=file)
        return self.run()

    def run(self):
        " load a chunk of text to be parsed and possibly executed"
        ret = None
        print '---- This is run! ' 
        while len(self.text)>0:
            s,nline,srcfile = self.text.pop()
#             print '   len(text) :=  %i (=%i+%i) ' % (len(self.text),
#                                                   len(self.text.front),
#                                                   len(self.text.back) )
#             print "next line = <<%s>> " % s
            
            if s == self.EOF:  return None

            if len(s)<=0: continue
            s.strip()
            if s.startswith('#'): continue
            print 'run-> eval <<%s>> ' % s
            ret = self.eval(s)
        return ret

    def eval(self,s):
        "evaluate a single tdl statement"
        ret = None
        token,code = self.compile(s)
        print 'compiled to ', token, code
        if token is None:
            for stmt in code:
                ret = self.expr_interp(stmt)
                if ret is not None and isinstance(ret,tuple):
                    key = ret[0]
                    if key == opcodes._print:
                        dest = ret[1] or sys.stdout
                        for i in ret[3]:   dest.write(i)
                        if ret[2]: dest.write("\n")
                    elif key == opcodes._continue:
                        print '? continue ' , ret

        print 'evaled to ', ret
        return ret
    
#        if x is not None:
#             ret = self.interpret(x, text=s)
#             self._status = True
#         return ret

    def get_input(self):
        return raw_input(self.prompt)

    def get_next_textline(self):

        filename = 'stdin'
        print 'This is get_next_textline!! ' 
        try:

            line,nline,filename  = self.text.pop()
        except IndexError:
            if self.interactive:
                print 'get more input: ', line, nline, filename
                line = self.get_input()
                nline = self.nline
            else:
                return (self.EOF,-1,'')
        self.nline = self.nline + 1
        return line,nline,filename

    def get_next_statement(self,s = None):
        " get and pre-process next program statement"
        if s is None:
            s,nline,fname = self.get_next_textline()
        s.strip()

        # print "GET NEXT STATEMENT ", s

        if s.startswith('#'): return ('','')
        # handle the case of triple quotes:
        #    join strings with a '\n' instead of ' '
        if s.find("'''")>-1   and not self.triple_dquote:
            self.triple_squote  = not self.triple_squote
        elif s.find('"""')>-1 and not self.triple_squote:
            self.triple_dquote  = not self.triple_dquote
        join = ' '
        if self.triple_dquote or self.triple_squote: join = '\n'

        if self.line_buff:  s = "%s%s%s" % (self.line_buff,join,s)

        jcom =find_unquoted_char(s,char='#')
        s = s[:jcom]

        n_parens = None
        # print '====Search for Closed Parens ', s
        # print '====   ', parens_matched(s)
        while n_parens != 0:
            n_parens = parens_matched(s)
            if n_parens < 0:  raise EvalError, 'syntax error: parens not matched: %s' % s
            if n_parens > 0:
                self.line_buff = s
                # print 'n_paremn>0?? ', n_parens, s
                s,nline,fname = self.get_next_textline()
                if s == self.EOF and not self.interactive:
                    break
                s = "%s%s%s" % (self.line_buff,join,s)
        #
        self.line_buff = ''
        self.triple_dquote = False
        self.triple_squote = False
        s = s.strip()
        if len(s)>1:
            # remove end-of-line comments
            jcom =find_unquoted_char(s,char='#')
            s = s[:jcom]
            # look for unquoted semi-colons (multiple statements per line)
            jsemi =find_unquoted_char(s,char=';')
            if jsemi<jcom:
                self.text.append((s[jsemi+1:],self.nline,self.infile))
                s = s[:jsemi]
            s.strip()
            # look for certain "keyword(x)" constructs:
            for j in ('if','elif', 'while', 'return', 'print'):
                if s.startswith("%s("% j): s = "%s %s" % (j,s[len(j):])
            # check for 'else:'
            if s in ('else:','try:','except:'): s = '%s :' % s[:-1]
        if len(s) < 1: s = ''
        # get first word:
        w = s.split()
        key = ''
        if len(w)>0: key = w[0].lower()
        return (s,key)

    def compile(self, expr= None):
        """main compilation of tdl statements: converts """
        print 'compile -> get next statement ', expr, len(self.text)
        s,key = self.get_next_statement(s=expr)
        print "compile get_next_st returned :key <%s> / statement= <%s> " % ( key, s)

        if s in ('','#'): return None
        if s.startswith('#'): return None
        if s == self.EOF: return self.EOF, []

        if key not in self.reserved_words:
            return None, [self.expr_compile(s)]

        elif key in ('del','print', 'return','break', 'continue'):
            # keywords that just get handled as special statements
            return key, [self.expr_compile(s)]

        # these keywords can never legally occur at a top-level compilation
        # and will always be found while the appropriate (if,while,for,def)
        # block is being processed
        if key in ('else','elif','endif','endfor','endwhile','enddef','endtry'):
            raise EvalError, 'syntax error: %s' % key

        # handle if statements, including if/elif/else/endif blocks
        elif key == 'if':
            print 'IF not implemented'
        elif key == 'try':
            print 'TRY not implemented'
        # def, for, and while blocks
        elif key in ('def', 'for', 'while'):
            print 'DEF/FOR/WHILE (%s) not implemented' % key


        # regular assignment / eval statement
        else:
            # check if command-like interpretation is reasonable
            # print 'check for command / assignment '
            print 'REGULAR assignment ', s
            try:
                next_char = s[len(key):].strip()[0]
            except:
                next_char = None
            if (key.find('(')==-1  and key.find(',')==-1 and key.find('=')==-1 and
                next_char != '='):
                if self.symbolTable.hasSymbol(key):
                    t = Command2Expr(s,symtable=self.symbolTable)
                    stack= self.expr_compile(t)
                    if stack[0] != opcodes.function:
                        raise EvalError, 'syntax error: weird error with commad '
                    stack[0] = opcodes.command
                    return  ['eval', stack, s]
            # wasn't a command!
            status,s1,s2 = split_delim(s,delim='=')
            if status == -1: return None
            if s2 == '': # Eval
                return  ['eval', self.expr_compile(s1), s]
            else: # Assignment
                stack = self.expr_compile(s1)
                tok = stack.pop(0)
                if tok not in (opcodes.variable, opcodes.array):
                    raise EvalError, 'syntax error: invalid assignment statement'
                if tok ==  opcodes.variable:  stack.insert(0,0)
                stack.insert(0,opcodes.symbol)
                return ['assign', stack,  self.expr_compile(s2), s]

        return ret



        ################
        
    def xxcompile(self, expr= None):
        """main compilation of tdl statements: converts a
        statement (or line of code) to executable tdl code.

        this relies heavily on Expression.compile to compile
        tdl (err, Python) statements to Python's AST format.

        The job here is to look for reserved keywords for control-flow
        code,  and build up a list of compiled (ie, AST-formatted)
        statements representing the blocks in the control flow.
        
        a tuple is returned of
             keyword, list of AST statements
        """
        print "compiling <%s>" %  expr
        s,key = self.get_next_statement(s=expr)
        print " :key <%s>" % key
        print " :s  <%s>"  % s
        if s in ('','#'): return None
        if s.startswith('#'): return None
        if s == self.EOF: return self.EOF, []

        if key not in self.reserved_words:
            return None, [self.expr_compile(s)]
        

        # these keywords can never legally occur at a top-level compilation
        # and will always be found while the appropriate (if,while,for,def)
        # block is being processed
        if key in ('else','elif','endif','endfor','endwhile','enddef','endtry'):
            raise EvalError, 'syntax error: %s' % key

        # handle if statements, including if/elif/else/endif blocks
        elif key == 'if':
            t = s[len(key):]
            status,s1,s2 = split_delim(t, delim=':')
            if status == -1:   return None
            ret = ['if']
            blockhead = self.expr_compile(s1)
            # block to execute
            if s2:   #  simple 'if x : y = x' type
                ret.append( (blockhead, [self.compile(s = s2)]))
            else:
                end = "endif"
                t = None
                block = []
                tmp   = []
                cond = [blockhead]
                else_seen = False
                while True:
                    sn,nextkey = self.get_next_statement()
                    if sn == '': continue
                    if sn == self.EOF:
                        raise EvalError, "end of file in 'if' block, file ='%s'" % self.infile 
                    if nextkey == 'endif':
                        block.append(tmp)
                        break
                    elif nextkey == 'elif':
                        if else_seen:
                            raise EvalError, 'syntax error: elif after else'
                        t = sn[len(nextkey):]
                        status,s1,s2 = split_delim(t, delim=':')
                        # print 'ELIF ', status , s1, s2
                        if status == -1:   return None
                        cond.append(self.expr_compile(s1))
                        block.append(tmp)
                        tmp = []
                    elif nextkey == 'else':
                        else_seen = True
                        cond.append(self.expr_compile('1'))
                        block.append(tmp)
                        tmp = []
                    else:
                        t = self.compile(s=sn)
                        if t is not None:  tmp.append(t)

                for i in zip(cond,block):  ret.append(i)
        # try: except: endtry:
        elif key == 'try':
            t = s[len(key):]
            status,s1,s2 = split_delim(t, delim=':')
            ret = ['try']
            blockhead = self.expr_compile('1')
            if s2:
                raise EvalError, 'syntax error: invalid try statement.'
            else:
                end = "endtry"
                t = None
                block = []
                tmp   = []
                cond = [blockhead]
                else_seen = False
                while True:
                    sn,nextkey = self.get_next_statement()
                    if sn == '': continue
                    if sn == self.EOF:
                        raise EvalError, "end of file in 'try' block, file ='%s'" % self.infile 
                    if nextkey == 'endtry':
                        block.append(tmp)
                        break
                    elif nextkey == 'except':
                        else_seen = True
                        cond.append(self.expr_compile('1'))
                        block.append(tmp)
                        tmp = []
                    else:
                        t = self.compile(s=sn)
                        if t is not None:  tmp.append(t)
                for i in zip(cond,block):  ret.append(i)
        #
        # def, for, and while blocks
        elif key in ('def', 'for', 'while'):
            t = s[len(key):]
            # print 'This is a ', key, ' t= ', t
            status,s1,s2 = split_delim(t, delim=':')
            if status == -1: return None
            ret = [key]
            # while blockhead is the simplest of all
            if key == 'while':
                blockhead = self.expr_compile(s1)
            elif key  == 'for':
                status,r1,r2 = split_delim(s1.replace(' in ',' @ '), delim='@')
                if status == -1: return None
                j = self.expr_compile(r1, reverse=True)
                x = j.pop()
                if len(j)!=1 or x != opcodes.variable:
                    raise EvalError, 'syntax error: invalid for A statement'
                blockhead = (j.pop(), self.expr_compile(r2))

            elif key == 'def':
                # there are three forms of def:
                #   def x = expr
                # and
                #   def x(args): return expr
                # and
                #   def x(args):
                #      block
                #   enddef
                # if status (holding position of ':') is < 4,
                # it cannot be either of the last two forms

                if status <4:
                    status,s1,s2 = split_delim(t, delim='=')
                    if status<1:
                        raise EvalError, 'syntax error: invalid def statement'
                    # it is of the simple form:   def x = expr
                    stack = self.expr_compile(s1, reverse=True)
                    tok = stack.pop()
                    if tok != opcodes.variable:
                        raise EvalError, 'syntax error: invalid def statement'
                    stack.append(opcodes.symbol)

                    return ['defvar',stack, self.expr_compile(s2),s2]
                else:
                    n1,n2 =  s1.find('('), s1.rfind(')')
                    if n1 <1 or n2<n1:
                        raise EvalError, 'syntax error: invalid def statement'

                    try:
                        t = self.expr_compile(s1,  reverse=True)
                    except:
                        raise EvalError, 'syntax error: invalid def statement'

                    ftype = t.pop() ;nargs = t.pop() ;  fname = t.pop()
                    if ftype != opcodes.function:
                        raise EvalError, 'syntax error: invalid def statement'

                    # construct tuple and keywords of function arguments
                    vargs = [] ; kws = {} ; iargs = 0; eq_seen = False
                    if n2>n1+1:
                        for i in split_list(s1[n1+1:n2]):
                            iargs = iargs+1
                            ieq = i.find('=')
                            if ieq == -1:
                                if eq_seen:
                                    raise EvalError, 'syntax error: invalid def statement 2'
                                j = self.expr_compile(i,reverse=True)
                                if j is not None and len(j)>1:
                                    x = j.pop()
                                    if len(j)!=1 or x != opcodes.variable:
                                        raise EvalError, 'syntax error: invalid def statement 1'
                                    vargs.append(j.pop())
                                else:
                                    iargs = iargs-1
                            else:
                                eq_seen = True
                                j = self.expr_compile(i[:ieq], reverse=True)
                                x = j.pop()
                                if len(j)!=1 or x != opcodes.variable:
                                    raise EvalError, 'syntax error: invalid def statement 1'
                                k = j.pop()
                                v = self.expr_eval(self.expr_compile(i[ieq+1:]))
                                kws[k] = v
                    if iargs != int(nargs):
                        raise EvalError, 'syntax error: invalid def statement 5'
                blockhead = [fname, tuple(vargs), kws]
            ret.append(blockhead)
            # block to execute
            if s2:   #  simple 'for i in arange(10): print i' type
                ret.append([self.compile(s = s2)])
            else:
                end = "end%s" % key
                t = 0 ; tmp = []
                while True:
                    sn,nextkey = self.get_next_statement()
                    if sn == '': continue
                    if sn == self.EOF:
                        raise EvalError, "end of file in '%s' block, file ='%s'" % (key,self.infile )
                    
                    if nextkey == end:
                        break
                    else:
                        t = self.compile(s=sn)
                        if t is not None:  tmp.append(t)
                ret.append(tmp)

        elif key in ('del','print', 'return'): # keywords that take a list
            s = s[len(key):].strip()
            for pars in (('(',')'),('[',']')):
                if (s.startswith(pars[0]) and s.endswith(pars[1])):
                    i = s.find(pars[1])
                    if i  == len(s)-1: s= s[1:len(s)-1]
            x = None
            ret.append(key.lower())
            if len(s)>0:
                ret.append( self.expr_compile("[%s]" % s))
            else:
                ret.append('')

        elif key in ('break', 'continue'):
            # break and continure are closely related,
            # and at parse stage do nothing.
            s = s[len(key):].strip()
            if len(s)>0:
                raise EvalError, 'syntax error: invalid %s statement' % key
            ret.append(key)

        # regular assignment / eval statement
        else:
            # check if command-like interpretation is reasonable
            # print 'check for command / assignment '
            # print 'REGULAR assignment ', s
            try:
                next_char = s[len(key):].strip()[0]
            except:
                next_char = None
            if (key.find('(')==-1  and key.find(',')==-1 and key.find('=')==-1 and
                next_char != '='):
                if self.symbolTable.hasFunction(key):
                    t = Command2Expr(s,symtable=self.symbolTable)
                    stack= self.expr_compile(t)
                    if stack[0] != opcodes.function:
                        raise EvalError, 'syntax error: weird error with commad '
                    stack[0] = opcodes.command
                    return  ['eval', stack, s]
            # wasn't a command!
            status,s1,s2 = split_delim(s,delim='=')
            if status == -1: return None
            if s2 == '': # Eval
                return  ['eval', self.expr_compile(s1), s]
            else: # Assignment
                stack = self.expr_compile(s1)
                tok = stack.pop(0)
                if tok not in (opcodes.variable, opcodes.array):
                    raise EvalError, 'syntax error: invalid assignment statement'
                if tok ==  opcodes.variable:  stack.insert(0,0)
                stack.insert(0,opcodes.symbol)
                return ['assign', stack,  self.expr_compile(s2), s]

        return ret

    def do_assign(self,left_hs,rhs):
        """ handle assignment statements """
        lhs  = left_hs[:]
        # the symbol name for the lhs will be looked for only
        # in the 'current group' unless fully qualified
        
        ndim_lhs = lhs.pop()
        varname  = lhs.pop()

        sym  = self.symbolTable.getSymbolLocalGroup(varname)
       
        # if the rhs evaluates to a symbol group, figure out where to place it
        # print 'Do Assign ', varname, sym, lhs, rhs

        if isGroup(rhs):
            return self.symbolTable.placeGroup(rhs,sym.name)

        if sym is None:
            self.raise_error('Cannot make assignment to %s??' % varname)
        for i in range(len(lhs)):
            if type(lhs[i]) == types.FloatType: lhs[i] = int(lhs[i])

        if isSymbol(sym):
            if sym.constant:
                self.raise_error('Cannot reassign value of %s' % sym.name)
            if len(lhs)>0 and  sym.type==symTypes.defvar:
                self.raise_error('Cannot assign to part of defined variable %s.' % varname)

        try:
            x = sym.value             # x = copy.deepcopy(sym.value)
        except:
            self.raise_error('Cannot make assignment %s' % sym.name)

        if len(lhs)==0:
            x = rhs
        elif len(lhs)==1:
            x[lhs[0]] = rhs
        else:
            x[tuple(lhs)] = rhs

        # 
        if sym.type in (symTypes.defvar,symTypes.defpro,symTypes.variable):
            sym.value  = x
            sym.code   = None
            sym.type   = symTypes.variable
            sym.constant = False

    def interpret(self,s,code=None,text=''):
        "interpret parsed code from compile"
        print 'Interpret: token= ', tok
        print 'Interpret: code ', code
        tok = s[0]
        if text != '': self.Expression.text=text

        # most common case: simple statement (non control-flow code)
        if tok is None:      return self.expr_interp(s[1])

        if tok == self.EOF:  return None
        try:
            tok=tok.lower()
        except AttributeError:
            pass

        if tok in self.__interrupts:
            self.interrupt = self.__interrupts.index(tok)
            if tok == 'return':
                self.retval = self.expr_eval(s[1])
                return
        elif tok == 'del':
            try:
                s[1].reverse()
                xtok = s[1].pop()
                nvar = s[1].pop()
            except:
                xtok = None
                
            if xtok != opcodes.list:
                raise EvalError, 'Invalid "del" statement'
            try:
                for i in range(nvar):
                    xtok = s[1].pop()
                    vname = s[1].pop()
                    if xtok != opcodes.variable:
                        raise EvalError, 'cannot delete %s ' % vname
                    self.symbolTable.delSymbol(vname)
            except:
                raise EvalError, 'Invalid "del" statement'
        elif tok == 'print':
            if len(s[1])>0:
                for i in self.expr_eval(s[1]):
                    self.output.write('%s ' % str(i))
            self.output.write('\n')
        elif tok == 'defvar':
            if len(s[1]) != 2 or s[1][1] != opcodes.symbol:
                raise EvalError, 'Invalid "def" statement'
            x = self.symbolTable.setDefVariable(s[1][0], s[2], s[3])
        elif tok == 'assign':
            lhs = self.expr_eval(s[1])
            rhs = self.expr_eval(s[2])
            self.do_assign(lhs,rhs)
        elif tok == 'eval':
            return self.expr_eval(s[1])
        elif tok == 'try':
            do_except = False
            for sx in  s[1][1]:
                try:
                    ret = self.interpret(sx)
                except:
                    do_except = True
                    break
            if do_except:
                for sx in s[2][1]:
                    ret = self.interpret(sx)

        elif tok == 'if':
            self.interrupt = 0
            found = False
            ret = None
            for i in s[1:] :
                cond,block = i[0],i[1]
                if self.expr_eval(cond):
                    for sx in block:
                        ret = self.interpret(sx,text)
                        if self.interrupt>0: break
                    if self.interrupt>1:  break
                    # self.interrupt = 0
                    return ret
        elif tok == 'while':
            self.interrupt = 0
            cond = self.expr_eval(s[1])
            while cond:
                for sx in s[2]:
                    ret = self.interpret(sx,text)
                    if self.interrupt>0: break
                if self.interrupt>1:  break
                self.interrupt = 0
                cond = self.expr_eval(s[1])
        elif tok == 'for':
            self.interrupt = 0
            iter_var = s[1][0]
            rhs = self.expr_eval(s[1][1])
            for x in rhs:
                self.symbolTable.setVariable(iter_var,x)
                for sx in s[2]:
                    ret = self.interpret(sx,text)
                    if self.interrupt>0: break
                if self.interrupt>1:  break
                self.interrupt = 0
        elif tok == 'def':
            #  look for docstring
            desc = None
            code = s[2]
            try:
                tx   = s[2][0]
            except IndexError:
                tx = ['']

            if (tx[0] == 'eval' and  len(tx[1])==2 and tx[1][0] == opcodes.string):
                desc = trimstring(tx[1][1])
                code = s[2][1:]

            self.symbolTable.setProcedure(s[1][0], code, desc=desc,
                                          args=s[1][1],kws=s[1][2])
        elif tok == self.EOF:
            return None
        else:
            raise EvalError, 'unknown evaluation error'

    ######

    def run_procedure(self,proc,args=None,kws=None):
        " run a user-created tdl procedure "
        if proc.type != symTypes.defpro:
            raise EvalError, 'invalid procedure'

        name = proc.name
        if len(args) != len(proc.args):
            raise EvalError, 'not enough arguments for procedure %s' % name

        name = "%s" % name.replace('.','@') # mangle module/function name 
        locgroup = self.symbolTable.LocalGroup
        modgroup = self.symbolTable.ModuleGroup

        group = self.symbolTable.addTempGroup(prefix=name,nlen=4)

        self.symbolTable.LocalGroup  = group.name
        self.symbolTable.ModuleGroup = proc.mod

        if group is None:
            raise EvalError, 'cannot run procedure %s (cannot create group??)' % name

        for k,v in zip(proc.args,args):
            self.symbolTable.setSymbol("%s.%s" % (group.name,k),v)

        kvals = {}
        kvals.update(proc.kws)
        kvals.update(kws)
        for k,v in kvals.items():
            self.symbolTable.setSymbol("%s.%s" % (group.name,k),v)

        ret = None
        try:
            for i in proc.code:
                self.interpret(i)
                if self.interrupt == 3:
                    ret = self.retval
                    self.interrupt = 0
                    break
        except:
            s = 'Error in procedure %s\n    %s' % (proc.name, i[-1])
            self.ShowError(msg=s,showtraceback=False)

            print proc.code
            
        try:
            if len(ret) == 1: ret= ret[0]
        except TypeError:
            pass

        self.symbolTable.delGroup(group.name)
        self.symbolTable.LocalGroup=locgroup
        self.symbolTable.ModuleGroup= modgroup
        return ret

    def ShowError(self,msg=None,showtraceback=True):
        " print error on exceptions"
        w = self.output.write
        if True: # try:
            w('************************************\n')
            w("%s\n" % msg)
            exctype, val, tb  = sys.exc_info()
            if not showtraceback: tb = None
            sys.excepthook(exctype,val,tb)
        else: # except:
            w('***Error printing exception error***\n')
            
