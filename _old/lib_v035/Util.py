# M. Newville Univ of Chicago (2006)
#
# -------------
# Modifications
# -------------
# 6-10-06 T2:
# - redid file_open as a class, can globally override open
#
# 4-30-06 MN:
# - improved list2array, show_list
#
# 4-29-06 T2:
# - added a file_open fcn
#
# 4-16-2006 T2:
# - added set_path and associated functions
#
# 1-3-06 T2:
#  - added mod_importer and a few other functions
#  - moved selftest into new TdlTest module.
#
##########################################################################
from Num import Num
import types
import exceptions
import sys
import string
import os

class ParseException(exceptions.Exception):
    def __init__(self,args=None):
        self.args = args

class EvalException(exceptions.Exception):
    def __init__(self,args=None):
        self.args = args
        
####
def verify_tdl(tdl,name='unknown',msg=''):
    if tdl is None:
        mout = "No tdl reference in python function '%s'" % (name)
        if msg != '': mout = "%s\n%s" % (mout,msg)
        raise RunTimeError, mout
        
def datalen(x):
    "return length of data for many datatypes"
    try:
        return len(x)
    except:
        return 1

def show_more(text,filename=None,writer=sys.stdout,pagesize=30,prefix=''):
    """ show lines of text in the style of more """
    txt = text[:]
    if type(txt)== types.StringType:  txt = txt.split('\n')        

    if len(txt) <1: return
    prompt = '== hit return for more, q to quit'
    ps = "%s (%%.2f%%%%) == " % prompt
    if filename: ps = "%s (%%.2f%%%%  of %s) == " % (prompt,filename)

    i = 0
    for i in range(len(txt)):
        if txt[i].endswith('\n'): writer.write("%s%s" % (prefix,txt[i]))
        else:                     writer.write("%s%s\n" % (prefix,txt[i]))
        i = i + 1
        if i % pagesize == 0:
            try:
                x = raw_input(ps %  (100.*i/len(txt)))
                if x in ('q','Q'): return
            except KeyboardInterrupt:
                writer.write("\n")
                return

def show_list(lst,ncol=None,textwidth=72):
    "formatted list of stuff in a list or dict"
    nmax = -1
    if len(lst) == 0: return ""
    if ncol is None:
        ncol = 1
        if type(lst) == types.ListType:        
            ncol = textwidth / (2 + max([len(i) for i in lst]))
    if ncol is None or ncol < 1: ncol = 1

    fmt = '%-'+str(int(textwidth/ncol))+'s  '

    if type(lst) == types.ListType:
        x = [str(c) for c in lst]
    elif type(lst) == types.DictType:
        x = [" %s= %s" % (k,str(lst[k])) for k in lst.keys()]
        
    pstr, num = '',1
    for c in x:
        t = fmt % str(c).rstrip()
        pstr = "%s%s" % (pstr, t)
        if num == ncol:
            pstr = "%s\n" % pstr
            num  = 0
        num = num + 1
    return  pstr

def find_unquoted_char(s,char='#'):
    "find character in a string, skipping over quoted text"
    i = 0 
    x = s.find(char)
    if x==-1: return len(s)
    use_triple = False
    while i < len(s):
        t = s[i]
        if t in ('"',"'"):
            if i < len(s)-3 and s[i:i+3] in ("'''",'"""'):  t = s[i:i+3]
            q = find_matching_quote(s[i:],quote=t)
            if q[0]: i = q[2]+i
        elif t==char:
            return i
        i = i + 1
    return i

def split_list(s,delim=','):
    "find character in a string, skipping over quoted text"
    i = 0
    quotes = {'"':'"',"'":"'","(":")","[":"]",'"""':'"""',"'''":"'''"}
    x = s.find(delim)
    if x==-1: return (s,)
    r = []
    i0 = 0
    while i < len(s):
        t = s[i]
        if t in quotes.keys():
            q = find_matching_quote(s[i:],quote=t,match=quotes[t])
            if q[0]: i = q[2]+i 
        elif t==delim:
            r.append(s[i0:i])
            i0 = i+1
        i = i + 1
    if i > i0: r.append(s[i0:])
    return r

def split_arg_str(s):
    "split s on ',' and on whitespace, keep key=val pairs together"
    args = []
    #first = s.split(',')
    first = split_list(s,',')
    for x in first:
        if x.find('=') == -1:
            #second = x.split()
            second = split_list(x,' ')
            for tok in second:
                y = tok.strip()
                if len(y) > 0:
                    args.append(y)
        else:
            y = x.strip()
            args.append(y)
    return args

def split_delim(s,delim='='):
    """ given a string of 'program text', split into parts around a
    single delimeter, such as an '=' sign for an assignment statement
    or a ':' for a for,if,def, or while statement.

    skips over matching quotes, and checks for matching brackets,
    braces,and parens.  will return status = -1 for incomplete statements.

    'bracket depth'=0, and also skips over matched single/double quotes.
    """
    opens = ['{','(','[','"',"'"]
    close = ['}',')',']','"',"'"]
    depth = [ 0, 0, 0, 0, 0]

    if not parens_matched(s): return (-1,s,'')
    i,idel,p,n,t = (0,0,None,None,None)
    while i < len(s):
        if i>1:        p = s[i-1]
        if i<len(s)-1: n = s[i+1]
        t = s[i]
        if t==delim and sum(depth)==0:
            if t!='=' or (t=='=' and n != '=' and idel == 0 and
                p not in ('!','>','<','=')):
                idel = i
                
        elif t in ('"',"'"):
            if i < len(s)-3 and s[i:i+3] in ("'''",'"""'):  t = s[i:i+3]
            q = find_matching_quote(s[i:],quote=t)
            if q[0]:
                i = q[2]+i
            else:
                j = opens.index(t)
                depth[j] = 1
        elif t in opens:
            j = opens.index(t)
            depth[j] = depth[j]+1
        elif t in close:
            j = close.index(t)
            depth[j] = depth[j]-1
        i = i + 1
    if sum(depth) != 0:
        return (-1,s,'')
    elif idel > 0 and idel<len(s):
        return (idel,s[:idel].strip(),s[idel+1:].strip())
    else:
        return (idel,s.strip(),'')

def unescape_string(s):
    escapes =(("\\n","\n"), ("\\r","\r"), ("\\a","\a"),("\\f","\f"),
              ("\\\\","\\"),
              ("\\b","\b"), ("\\v","\v"), ("\\t","\t"))
    for i,j in escapes: s = s.replace(i,j)
    return s

def trimstring(s,use_raw=False):
    " trim leading 'quotes' from string variables"
    if type(s) != types.StringType: return s
    if (s.startswith("'''")   and s.endswith("'''")):  return r'%s' % s[3:-3]
    elif (s.startswith('"""') and s.endswith('"""')):  return unescape_string(s[3:-3])
    elif (s.startswith("'")   and s.endswith("'")):    return r'%s' % s[1:-1]
    elif (s.startswith("\"")  and s.endswith("\"")):   return unescape_string(s[1:-1])
    return s


def _isnumericarray(x):
    """returns whether value (potentially nested list) can be coerced to numerical array.
    note that this insists on using numpy for numerical arrays, not record arrays.
    """
    for i in x:
        ok = True
        if type(i) in (types.ListType,types.TupleType):
            ok = _isnumericarray(i)
        elif type(i) not in (types.FloatType,types.ComplexType,types.IntType,types.LongType):
            ok = False
        if not ok: return False
    return True

def list2array(s):
    """ attempt to convert a list to a NumPy array.
    Returns original list if the conversion is not possible"""
    if type(s) == types.ListType:
        try:
            if _isnumericarray(s):
                s = Num.array(s)
        except (TypeError,AttributeError):
            pass
    return s

def strip_ending_comma(s):
    " strip ending comma from a string"
    if type(s)== types.StringType and s.endswith(','):  return s[:-1]
    return s

def find_matching_quote(s,quote='"',match=None):
    """
    find positions of matching quotes in a string
    """
    esc = "\\"
    j = s.find(quote)
    if match is None: match = quote
    if j > -1 and s[j:j+len(quote)] == quote:
        p,k = (None,j)
        while k < j+len(s[j+1:]):
            k = k+1
            if s[k:k+len(match)] == match and p != esc:
                return True,j,k+len(match)-1
            p = s[k:k+1]
    return False,j,len(s)

def parens_matched(s):
    """ given a string of 'program text',
    skips over matching quotes, and checks for matching brackets,
    braces,and parens.  returns True / False
    """
    delims = {'{':1,'(':1,'[':1, '}':-1,')':-1,']':-1}
    b = delims.keys()
    i = 0
    depth = 0
    while i < len(s):
        t = s[i]
        if t in ('"',"'"):
            if i < len(s)-3 and s[i:i+3] in ("'''",'"""'):  t = s[i:i+3]
            q = find_matching_quote(s[i:],quote=t)
            # print q
            if q[0]:
                i = q[2]+i

            else: # if matching quote not found, string is not complete
                return False
        elif t in b:
            depth = depth + delims[t]
        i = i + 1
    return (depth == 0)

def Command2Expr(s,symtable=None):
    """ convert command like syntax to function syntax:
    in a command, commas are optional / superfluous, and
    there may be extra space are key=val arguments:

       cmd arg1  arg2, key1= val key2= val
    =>
       cmd(arg1, arg2, key1=val, key2=val)
    """
    key = s.split()[0].lower()
    s = s[len(key):]
    words = split_list(s, delim = ' ')

    tmp = []
    for word in words:
        if word.startswith(','): word = word[1:]
        if word.endswith(','):   word = word[:-1]
        word.strip()
        if word != '':
            if word.find(',')==-1: 
                tmp.append(word)
            else:
                subwords = split_list(word,delim=',')
                for s in subwords: tmp.append(s)
    s = ' '.join(tmp)

    if s.find('=') > -1:
        s = '='.join([i.strip() for i in split_list(s, delim='=')])

    # if a symboltable is included, decide which strings need quotes
    if symtable is not None:
        def needs_quote(w):
            'word is not a literal string, number, or a named symbol'
            if  (w.startswith('"') or w.startswith("'")): return False
            try:
                x = float(w)
                return False
            except ValueError:
                pass
            return not symtable.hasSymbol(w) 

        words = split_list(s,delim=' ')
        tmp = []
        for i in words:
            out = i
            if i.find('=')>0:
                k,v = i.split('=')
                if needs_quote(v): out = '%s="%s"'% (k,v)
            elif needs_quote(i):
                out ='"%s"' % i
            tmp.append(out)
        s = ' '.join(tmp)
    #
    s = ', '.join([i.strip() for i in split_list(s, delim=' ')])    
    expr = "%s(%s)" % (key,s)
    return expr

def int2bin(x):
    """ convert integer to list of booleans: gauranteed to return 12 'bits'"""
    keys = ('000','001','010','011','100','101','110','111')
    o = []
    for i in list(''.join([keys[int(i)] for i in oct(x)])):
        o.append(i=='1')
    o.reverse()
    if len(o)<= 12:
        for i in range(12-len(o)): o.append(False)
    return o

def mod_import(name):
    """
    wrapper for imports/reloads
    given a module name, try reloading the module
    given a string name, try loading the given module
    Note if a str is passed the returned mod will point to
    the mod at the end of a mod list, eg:
    mod = mod_import('a.b.c')
    mod -> c
    """
    if not name:
        print "No module name"
        return(None)
    if type(name) == types.ModuleType:
        try:
            reload(name)
            return name
        except ImportError:
            s = 'Error loading module %s' % name
            PrintExceptErr(s) 
            return None
    elif type(name) == types.StringType:
        try: 
            mod = __import__(name)
            components = name.split('.')
            for comp in components[1:]:
                mod = getattr(mod, comp)
            return mod 
        except ImportError:
            s = 'Error loading module %s' % name
            PrintExceptErr(s) 
            return None
    else:
        return(name)

def PrintExceptErr(err_str,print_trace=True):
    " print error on exceptions" 
    try:
        print '\n***********************************'
        print err_str
        print 'Exception:', sys.exc_type
        xx, yy, zz = sys.exc_info()
        if print_trace == False: zz = ''
        sys.excepthook(xx,yy,zz)
        print '***********************************\n'
    except:
        print '*****Error printing exception error******'

def PrintShortExcept(err_str):
    " print error on exceptions" 
    try:
        print '\n***********************************'
        print err_str
        xx, yy, zz = sys.exc_info()
        sys.excepthook(xx,yy,None)
        print '***********************************\n'
    except:
        print '*****Error printing exception error******'

def set_path(pth=None,recurse=False,verbose=False,clean=True):
    "modify or return python path"
    if type(pth) == types.StringType: pth = pth.strip()
    if not pth:
        return 
    pth = os.path.abspath(pth)
    if os.path.exists(pth):
        if verbose: print 'add->', pth
        if pth not in sys.path: sys.path.append(pth)
        if recurse == True:
            dirs = sub_dirs(pth,skip_txt='.')
            for d in dirs:
                if d not in sys.path:
                    if os.path.exists(d):
                        if verbose: print 'add->', d
                        sys.path.append(d)
    else:
        if verbose: print "Path '%s' doesnt exist" % pth
    if clean: clean_path()
    return

def clean_path():
    temp = []
    for p in sys.path:
        if p not in temp:
            temp.append(p)
    # temp.sort()
    sys.path = temp

def sub_dirs(pth,skip_txt=None):
    sub_dirs = []
    if os.path.exists(pth):
        for root,dirs,files in os.walk(pth):
            for d in dirs:
                temp = os.path.abspath(os.path.join(root,d))
                if skip_txt:
                    if skip_txt not in temp:
                        sub_dirs.append(temp)
                else:
                    sub_dirs.append(temp)
    return sub_dirs

#def file_open(fname,default_path=None,**kw):
#    "open a file using default path if passed"
#    # two cases
#    # 1. fname has full path (or rel path to cwd), or file is in cwd -> dont use def_path
#    # 2. fname has rel path (or none) to default path -> join def_path and fname
#    print "hello file open: %s" % default_path
#    if os.path.isfile(fname):
#        #return open(fname)
#        return file(fname,**kw)
#    elif default_path:
#        fname = os.path.join(default_path,fname)
#        #return open(fname)
#        return file(fname,**kw)
#    else:
#        raise IOError, "Could not open file '%s'" % (fname)

class file_open:
    """open a file using default path.  The default path may be
    passed as a kw or determined from a symbol table if present.
    This class is intended to override the builtin function open
    The equivalent operation file is used here to return a file
    object.  Note do not override the builtin function 'file' with this class
    or you'll enter an infinite loop!!!  
     * two cases for finding the file:
     1. fname has full path (or rel path to cwd), or file is in cwd
        -> dont use def_path
     2. fname has rel path (or none) to default path
        -> join default_path and fname
    """
    def __init__(self,sym=None,file_path="_sys.work"):
        self.sym=sym
        self.file_path=file_path
    def __call__(self,*args,**kw):
        if kw.has_key("default_path"):
            default_path = kw.pop("default_path")
        elif self.sym and self.file_path:
            #default_path = self.sym.getSymbolValue("_sys.work")
            default_path = self.sym.getVariableValue(self.file_path)
        else:
            default_path = None
        #print default_path
        #print kw
        #print args
        if len(args) > 0:
            fname = args[0]
            args = args[1:]
        else:
            raise IOError, "No file name given"

        if len(args) > 0:
            if args[0] == 'w':
                return file(fname,*args,**kw) 
        
        if os.path.isfile(fname):
            #return open(fname)
            return file(fname,*args,**kw)
        elif default_path:
            fname = os.path.join(default_path,fname)
            #return open(fname)
            return file(fname,*args,**kw)
        else:
            raise IOError, "Could not open file '%s'" % (fname)

###################################################################
def testCommand2Expr():
    import Symbol
    s  = Symbol.SymbolTable()
    s.addVariable('avar',3)
    s.addVariable('bvar',3)
    s.addVariable('x',Num.arange(10))
    print Command2Expr('func avar, bvar cvar, dvar=3',symtable=s)
    print Command2Expr('func *.*  cvar, dvar=3',symtable=s)
    print Command2Expr('plot x x*2',symtable=s)
    

if __name__ == '__main__':
    print 'tdl utility functions.'
    
