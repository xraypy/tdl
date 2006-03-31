# M. Newville Univ of Chicago (2006)
#
# -------------
# Modifications
# -------------
#
# * 1-3-06 T2:
#  - added mod_importer and a few other functions
#  - moved selftest into new TdlTest module.
#
##########################################################################
from Num import Num
import types
import exceptions
import sys
import string

class ParseException(exceptions.Exception):
    def __init__(self,args=None):
        self.args = args

class EvalException(exceptions.Exception):
    def __init__(self,args=None):
        self.args = args
        
##
##
def datalen(x):
    "return length of data for many datatypes"
    try:
        return len(x)
    except:
        return 1



def show_more(text,filename=None,writer=sys.stdout,pagesize=24):
    """ show lines of text in the style of more """
    txt = text[:]
    if type(txt)== types.StringType:  txt = txt.split('\n')        

    if len(txt) <1: return
    prompt = '== hit return for more, q to quit'
    ps = "%s (%%.2f%%%%) == " % prompt
    if filename: ps = "%s (%%.2f%%%%  of %s) == " % (prompt,filename)
         
    i = 0
    for i in range(len(txt)):
        if txt[i].endswith('\n'): writer.write(txt[i])
        else:                     writer.write("%s\n" % txt[i])
        i = i + 1
        if i % pagesize == 0:
            try:
                x = raw_input(ps %  (100.*i/len(txt)))
                if x in ('q','Q'): return
            except KeyboardInterrupt:
                writer.write("\n")
                return

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
              ("\\b","\b"), ("\\v","\v"), ("\\t","\t"))
    for i,j in escapes: s = s.replace(i,j)
    return s

def trimstring(s):
    " trim leading 'quotes' from string variables"
    # print 'trim string inp: ', s

    if type(s) != types.StringType: return s
    if  ((s.startswith("'''")  and s.endswith("'''")) or
         (s.startswith('"""')  and s.endswith('"""'))):
        s = s[3:-3]
    elif ((s.startswith("'")  and s.endswith("'")) or
          (s.startswith("\"") and s.endswith("\"")) ):
        s = s[1:-1]
    return unescape_string(s)


def list2array(s):
    """ attempt to convert a list to a NumPy array.
    Returns original list if the conversion is not possible"""
    if type(s) == types.ListType:
        try:
            for i in s:  x = abs(i)
            s = Num.array(s)
        except TypeError:
            pass
    return s

def strip_ending_comma(s):
    " strip ending comma from a string"
    if type(s)==types.StringType and s.endswith(','):  return s[:-1]
    return s

def find_matching_quote(s,quote='"',match=None):
    """
    find positions of matching quotes in a string
    """
    esc = "\\"
    j = s.find(quote)
    if match == None: match = quote
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

def Command2Expr(key, s):
    """ convert command like syntax to function syntax:
    in a command, commas are optional / superfluous, and
    there may be extra space are key=val arguments:

       cmd arg1  arg2, key1= val key2= val
    =>
       cmd(arg1, arg2, key1=val, key2=val)
    """
    s = s[len(key):]
    a = split_list(s, delim = ' ')
    for i in range(len(a)):
        if a[i].startswith(','): a[i] = a[i][1:]
        if a[i].endswith(','):   a[i] = a[i][:-1]
        a[i].strip()
    b = []
    for i in a:
        if i.strip() != '': b.append(i)
    s = ' '.join(b)
    if s.find('=') > -1:
        s = '='.join([i.strip() for i in split_list(s, delim='=')])

    s = ', '.join([i.strip() for i in split_list(s, delim=' ')])    
    return "%s(%s)" % (key,s)
    
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
        except:
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
        except:
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


def show_list(lst,ncol=None,textwidth=80):
    "formatted list of stuff in a list or dict"
    nmax = -1
    if len(lst) == 0: return ""
    if ncol == None:
        ncol = 1
        if type(lst) == types.ListType:        
            ncol = textwidth / (2 + max([len(i) for i in lst]))

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



if __name__ == '__main__':
    print 'tdl utility functions.'
