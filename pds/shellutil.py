"""
Utilities for shell programs

Authors/Modifications:
---------------------
* M. Newville (2006)
* Modified from tdl-revision 226 for use with pds
  shell program Tom Trainor (tptrainor@alaska.edu)

"""
##########################################################################

import types
import time
import exceptions
import sys
import string
import os
import copy

##########################################################################
#  Data Utilities
##########################################################################
class Group:
    """Generic group"""
    def __init__(self,):
        pass
    def __repr__(self,):
        attr = []
        for a in dir(self):
            if a[0].isalpha():
                attr.append(a)
        if len(attr) > 0:
           lout = show_list(attr,textwidth=82)
        else:
            lout = 'Group is empty'
        return lout

def is_group(g):
    return isinstance(g,Group)

##########################################################################
def datalen(x):
    """
    return length of data for many datatypes
    """
    try:
        return len(x)
    except:
        return 1

##########################################################################
def _isnumericarray(x):
    """
    returns whether value (potentially nested list) can be coerced to
    numerical array. note that this insists on using numpy for numerical
    arrays, not record arrays.
    """
    for i in x:
        ok = True
        if type(i) in (types.ListType,types.TupleType):
            ok = _isnumericarray(i)
        elif type(i) not in (types.FloatType,types.ComplexType,
                             types.IntType,types.LongType):
            ok = False
        if not ok: return False
    return True

##########################################################################
def list2array(s):
    """
    Attempt to convert a list to a numpy array.
    Returns original list if the conversion is not possible
    """
    import numpy as num
    if type(s) == types.ListType:
        try:
            if _isnumericarray(s): return num.array(s)
        except:
            pass
    return s

##########################################################################
def str2list(s,conv=float):
    """
    convert a string to a list
    """
    if s == None: return []
    if s == 'None': return []
    if type(s) == types.StringType:
        s = s.strip()
        if s[0] != '[': s = '[' + s
        if s[-1] != ']': s = s + ']'
    s = eval(s)
    if len(s) == 0: return []

    #return map(conv, s)
    v = []
    for l in s:
        if conv == float:
            n = float(l)
            v.append(n)
        elif conv == int:
            n = int(l)
            v.append(n)
    return v

##########################################################################
def int2bin(x):
    """
    convert integer to list of booleans: gauranteed to return 12 'bits'
    """
    keys = ('000','001','010','011','100','101','110','111')
    o = []
    for i in list(''.join([keys[int(i)] for i in oct(x)])):
        o.append(i=='1')
    o.reverse()
    if len(o)<= 12:
        for i in range(12-len(o)): o.append(False)
    return o

##########################################################################

##########################################################################
#  Text/String/Parsing Utilities
##########################################################################

##########################################################################
def show_more(text,filename=None,writer=None,pagesize=30,prefix=''):
    """
    show lines of text in the style of more
    """
    txt = text[:]
    if type(txt)== types.StringType:  txt = txt.split('\n')

    if len(txt) <1: return
    prompt = '== hit return for more, q to quit'
    ps = "%s (%%.2f%%%%) == " % prompt
    if filename: ps = "%s (%%.2f%%%%  of %s) == " % (prompt,filename)

    if writer == None:
        writer = sys.stdout

    i = 0
    for i in range(len(txt)):
        if txt[i].endswith('\n'):
            writer.write("%s%s" % (prefix,txt[i]))
        else:
            writer.write("%s%s\n" % (prefix,txt[i]))
        i = i + 1
        if i % pagesize == 0:
            try:
                x = raw_input(ps %  (100.*i/len(txt)))
                if x in ('q','Q'): return
            except KeyboardInterrupt:
                writer.write("\n")
                return

##########################################################################
def show_list(lst,ncol=None,textwidth=72):
    """
    formatted list of stuff in a list or dict
    """
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

##########################################################################
def find_unquoted_char(s,char='#'):
    """
    find character in a string, skipping over quoted text
    """
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

##########################################################################
def split_list(s,delim=','):
    """
    split a string at delim, skipping over occurances of delim
    if they appear in quoted text
    """
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

##########################################################################
def split_cmd_line(s):
    """
    split a command line based on ';'
    """
    s = str(s)
    if s.find(';') == -1:
        return([s])
    tmp = s.split(';')
    lines = []
    for l in tmp: lines.append(l.strip())
    return lines

##########################################################################
def split_args(s):
    """
    split 's' on ',' and on whitespace, keep key=val pairs together
    """
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

##########################################################################
def split_delim(s,delim='='):
    """
    Given a string of 'program text', split into parts around a
    single delimeter, such as an '=' sign for an assignment statement
    or a ':' for a for,if,def, or while statement.

    skips over matching quotes, and checks for matching brackets,
    braces,and parens.  will return status = -1 for incomplete statements.

    'bracket depth'=0, and also skips over matched single/double quotes.
    """
    opens = ['{','(','[','"',"'"]
    close = ['}',')',']','"',"'"]
    depth = [ 0, 0, 0, 0, 0]

    if parens_matched(s) != 0: return (-1,s,'')
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
                print 'Split Delim: ', s, ' :: ', t
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

##########################################################################
def unescape_string(s):
    """
    replace escaped characters
    """
    escapes =(("\\n","\n"), ("\\r","\r"), ("\\a","\a"),("\\f","\f"),
              ("\\\\","\\"), ("\\b","\b"), ("\\v","\v"), ("\\t","\t"))
    for i,j in escapes: s = s.replace(i,j)
    return s

##########################################################################
def escape_string(s):
    """
    escape control characters
    """
    escapes =(("\n","\\n"), ("\r","\\r"), ("\a","\\a"),("\f","\\f"),
              ("\\","\\\\"), ("\b","\\b"), ("\v","\\v"), ("\t","\\t"))
    for i,j in escapes: s = s.replace(i,j)
    return s

##########################################################################
def trimstring(s):
    """
    trim leading 'quotes' from string variables
    """
    if type(s) != types.StringType: return s
    if (s.startswith("'''")   and s.endswith("'''")):  return r'%s' % s[3:-3]
    elif (s.startswith('"""') and s.endswith('"""')):  return unescape_string(s[3:-3])
    elif (s.startswith("'")   and s.endswith("'")):    return r'%s' % s[1:-1]
    elif (s.startswith("\"")  and s.endswith("\"")):   return unescape_string(s[1:-1])
    return s

##########################################################################
def strip_ending_comma(s):
    """
    strip ending comma from a string
    """
    if type(s)== types.StringType and s.endswith(','):  return s[:-1]
    return s

##########################################################################
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

##########################################################################
def parens_matched(s):
    """
    Given a string of 'program text',
    skips over matching quotes, and checks for matching brackets,
    braces,and parens.  returns:
         0  for completely matched string
        >0  for a string needing closing parens
        <0  for a string with extra closing parens (ie, a syntax error)
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
                return 1
        elif t in b:
            depth = depth + delims[t]
        i = i + 1
    return depth

##########################################################################
def command2expr(s,symtable=None):
    """
    Convert command like syntax to function syntax:
    in a command, commas are optional / superfluous, and
    there may be extra space in key=val arguments:

       cmd arg1  arg2, key1= val key2= val
    =>
       cmd(arg1, arg2, key1=val, key2=val)
    """
    key   = s.split()[0].lower()
    s     = s[len(key):]
    if len(s.strip()) == 0: return "%s()" % key

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
            """
            word is not a literal string, number, or a named symbol
            Note if string starts with a paren of sorts we assume its
            a raw data type.
            """
            if  (w.startswith('"') or w.startswith("'")):
                return False
            if  (w.startswith('[') or w.startswith("(") or w.startswith("{")):
                return False
            if  (w.startswith('.') ):
                return True
            try:
                x = float(w)
                return False
            except ValueError:
                pass
            return not symtable.has_symbol(w)
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

    # final expr
    s = ', '.join([i.strip() for i in split_list(s, delim=' ')])
    expr = "%s(%s)" % (key,s)

    return expr

##########################################################################

##########################################################################
#  Pickle Utilities
##########################################################################

##########################################################################
def pickle_1(data,fname,check=True):
    """
    write data to a pickle file
    data is assumed to be a dictionary!
    """
    import cPickle
    if data == None: return

    # see if data is pickleable
    if check:
        for dname in data.keys():
            # call save if has it
            if hasattr(data[dname],'__save__'):
                data[dname].__save__()
            try:
                # try to pickle to a string
                pstr = cPickle.dumps(data[dname])
            except:
                data.pop(dname)
                print "Warning object '%s' could not be pickled" % dname
        if len(data) == 0: return

    # pickle the data
    d = {'save_version': 1,
         'title':        'SaveSet',
         'pds_version':  1.0,
         'os_name':      os.name,
         'os_environ':   os.environ,
         'timestamp':    time.time(),
         'data':         data}
    f = open(fname,'w')
    cPickle.dump(d,f)
    f.close()

def unpickle_1(fname):
    """
    read data made by pickle 1
    """
    import cPickle

    if not os.path.isfile(fname):
        print "File error: cannot find file '%s' " % fname
        return None

    f = open(fname)
    d = cPickle.load(f)
    f.close()
    isOK = True
    try:
        vers = d['save_version']
        title= d['title']
        pver = d['pds_version']
        data = d['data']
        isOk = isOK and (d['title'] == 'SaveSet')
        isOk = isOK and (int(d['save_version']) >= 1.0)
    except:
        isOK = False

    if isOK and data != None:
        return data
    else:
        print "'%s' is not a proper save file" % fname
        return None

##########################################################################
def pickle_2(data,fname):
    """
    write data to a pickle file
    data is assumed to be a dictionary!
    this is untested!!
    """
    import cPickle
    if data == None: return

    # get pickleable data
    pdata = {}
    for dname in data.keys():
        # call save if has it
        if hasattr(data[dname],'__save__'):
            data[dname].__save__()
        try:
            # try to pickle to a string using
            # protocol 0
            pstr = cPickle.dumps(data[dname],0)
            pdata.update({dname:pstr})
        except:
            print "Warning object '%s' could not be pickled" % dname

    if len(pdata) == 0: return None
    f = open(fname,'w')
    for (dname,pstr) in pdata.items():
        s = "\n\n====>%s\n%s\n\n" % (dname,pstr)
        f.write(s)
    f.write("\n\n====>END\n")
    f.close()

def unpickle_2(fname):
    """
    read data made by pickle 2
    this is untested!!!
    """
    import cPickle

    if not os.path.isfile(fname):
        print 'file error: cannot find file %s ' % fname
        return None

    f = open(fname)
    lines = f.readlines()
    f.close()

    # get all the pickle strings
    pstrs = {}
    dname = ''
    pstr = ''
    for line in lines:
        if len(line.strip()) == 0:
            pass
        if line.startswith("====>"):
            if (len(dname)>0) and (len(pstr) > 0):
                pstrs.update({copy.copy(dname):copy.copy(pstr)})
            dname = line.split("====>")[1]
            pstr = ''
        else:
            pstr = pstr + line + '\n'
        if dname == 'END':
            break

    # pickle the data
    pdata = {}
    for (dname,pstr) in pstrs.items():
        try:
            dat = cPickle.loads(pstr)
            pdata.update({dname:dat})
        except:
            print "Could not restore object '%s'" % dname

    if len(pdata) > 0:
        return pdata
    else:
        return None

##########################################################################

##########################################################################
#  Module/Errors/Path/File Utilities
##########################################################################

##########################################################################
def mod_import(name,debug=False):
    """
    Wrapper for imports/reloads

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
            if debug:
                PrintExceptErr(s)
            else:
                print s
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
            if debug:
                PrintExceptErr(s)
            else:
                print s
            return None
    else:
        return(name)

##########################################################################
class  PrintExceptErr:
    """
    print error on exceptions
    """
    def __init__(self,err_str,print_trace=True):
        print '\n***********************************'
        #print 'PrintExceptErr', err_str
        try:
            print err_str
            print 'Error:', sys.exc_type
            xx, yy, zz = sys.exc_info()
            if print_trace == False: zz = ''
            sys.excepthook(xx,yy,zz)
        except:
            print '  Error printing exception error '
        print '***********************************\n'

##########################################################################
class  PrintShortExcept:
    """
    print error on exceptions
    """
    def __init__(self,err_str):
        print '\n***********************************'
        #print 'PrintShortExcept', err_str
        try:
            print err_str
            xx, yy, zz = sys.exc_info()
            sys.excepthook(xx,yy,None)
        except:
            print '  Error printing exception error '
        print '***********************************\n'

##########################################################################
def set_path(pth=None,recurse=False,verbose=False,clean=True):
    """
    modify or return python path
    """
    if pth is None:
        return
    if isinstance(pth, str):
        pth = pth.strip()
    if pth == ".":
        if '.' not in sys.path:
            sys.path.append('.')
    else:
        pth = os.path.abspath(pth)
        if os.path.exists(pth):
            if verbose: print '    add->', pth
            if pth not in sys.path: sys.path.append(pth)
            if recurse == True:
                dirs = sub_dirs(pth)
                for d in dirs:
                    if d not in sys.path:
                        if os.path.exists(d):
                            if verbose: print '    add->', d
                            sys.path.append(d)
        else:
            #if verbose: print "Path '%s' doesnt exist" % pth
            pass # print "Set path warning cant find path: '%s'" % pth
    if clean: clean_path()
    return

##########################################################################
def clean_path():
    """
    sort path strings
    """
    temp = []
    for p in sys.path:
        if p not in temp:
            temp.append(p)
    temp.sort()
    sys.path = temp

##########################################################################
def sub_dirs(pth,alp_only=True):
    """
    Find sub dirs
    """
    sub_dirs = []
    if os.path.exists(pth):
        for root,dirs,files in os.walk(pth):
            for d in dirs:
                skip = False
                p = os.path.abspath(os.path.join(root,d))
                if alp_only:
                    tmp = p.split(os.sep)
                    for t in tmp:
                        if len(t)>0 and not t[0].isalpha():
                            skip = True
                            break
                if not skip:
                    sub_dirs.append(p)
    return sub_dirs

##########################################################################
class FileOpen:
    """
    Open a file using default path.

    The default path may be passed as a kw or determined from a symbol
    table if present.

    This class can be used similiar to the builtin function open
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
            default_path = self.sym.get_symbol(self.file_path)
            #default_path = self.sym.getVariableValue(self.file_path)
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

        if os.path.isfile(fname):
            #return open(fname)
            return file(fname,*args,**kw)
        elif default_path:
            fname = os.path.join(default_path,fname)
            #return open(fname)
            return file(fname,*args,**kw)
        else:
            raise IOError, "Could not open file '%s'" % (fname)

##########################################################################

##########################################################################
#  Get Line / Menu Utilities
##########################################################################

##########################################################################
def gtline(pterm='>>',p0=None,p1=None,p2=None,valid=None,
           default=None,rettype=None,retry=True):
    """
    Fancy get line function

    Arguments:
    ----------
    * pterm   = prompt terminator
    * p0      = line to print prior to actual prompt
    * p1      = prompt
    * p2      = prompt on retry
    * valid   = list of valid entries
    * default = default ret value
    * rettype = typing function, ie function to force ret type
    * retry   = retry on failure?

    """
    if p0 != None:
        print p0
    if p1 != None:
        prompt = str(p1)
    else:
        prompt = ''
    ret = _gtline(prompt,default,pterm,rettype,valid)
    if (ret == None) and (retry == True):
        if p2 != None:
            prompt = str(p2)
        else:
            prompt = ''
        while 1:
            ret = _gtline(prompt,default,pterm,rettype,valid)
            if ret != None:
                break
    return (ret)

def _gtline(prompt='',default=None,pterm='>>',rettype=None,valid=None):
    if default != None:
        prompt = "%s (%s)%s" % (prompt,str(default),str(pterm))
    else:
        prompt = "%s%s" % (prompt,str(pterm))
    ret = raw_input(prompt)
    ret = ret.strip()
    if len(ret) > 0:
        if rettype != None:
            try:
                ret = rettype(ret)
            except:
                ret = None
    elif default != None:
        ret = default
    else:
        ret = None
    if valid != None:
        if ret not in valid:
            ret = None
    return (ret)

##########################################################################
def get_str(prompt=None,default=None,valid=None,nlprompt=True):
    """
    get a string from prompt
    """
    if prompt == None:
        prompt = ""
    if default != None:
        if valid == None:
            prompt = "%s (%s)" % (prompt,default)
        elif default in valid:
            prompt = "%s (%s)" % (prompt,default)
        else:
            default = None
    prompt = prompt + ">>"
    if nlprompt==True: prompt = prompt + "\n"
    while 1:
        try:
            ret = raw_input(prompt)
            ret = ret.strip()
            if len(ret) == 0 and default != None:
                return default
            elif valid != None:
                if ret in valid:
                    return ret
                else:
                    print 'Please enter a valid response'
                    print '-->', valid
            else:
                return ret
        except:
            pass

##########################################################################
def get_yn(prompt=None,default=None,nlprompt=True):
    """
    get a yes/no response from the prompt
    """
    ylbl=['yes','y','1']
    nlbl=['no','n','0']
    if prompt == None:
        prompt = "Enter Yes/No"
    if default != None:
        default = default.lower()
        if default in ylbl or default in nlbl:
            prompt = "%s (%s)" % (prompt,default)
        else:
            default = None
    prompt = prompt + ">>"
    if nlprompt==True: prompt = prompt + "\n"
    while 1:
        try:
            ret = raw_input(prompt)
            if ret.lower() in ylbl:
                return 'yes'
            elif ret.lower() in nlbl:
                return 'no'
            elif len(ret.strip()) == 0 and default != None:
                if default in ylbl:
                    return 'yes'
                elif default in nlbl:
                    return 'no'
            else:
                print 'Please enter Yes (%s) or No (%s)' % (str(ylbl),
                                                            str(nlbl))
        except:
            pass

##########################################################################
def get_tf(prompt=None,default=None,nlprompt=True):
    """
    get a true/false response from the prompt
    """
    tlbl=['true','t','1']
    flbl=['false','f','0']
    if prompt == None:
        prompt = "Enter True/False"
    if default != None:
        if type(default) == types.StringType:
            default = default.lower()
            if default in tlbl:
                default = True
            elif default in flbl:
                default = False
            else:
                default = None
        if default == True:
            prompt = "%s (%s)" % (prompt,'True')
        elif default == False:
            prompt = "%s (%s)" % (prompt,'False')
        else:
            default = None
    prompt = prompt + ">>"
    if nlprompt==True: prompt = prompt + "\n"
    while 1:
        try:
            ret = raw_input(prompt)
            if len(ret.strip()) == 0 and default != None:
                return default
            elif ret.lower() in tlbl:
                return True
            elif ret.lower() in flbl:
                return False
            else:
                print 'Enter True (%s) or False (%s)' % (str(tlbl),
                                                         str(flbl))
        except:
            pass

##########################################################################
def get_int(prompt=None,default=None,valid=[],min=None,max=None,nlprompt=True):
    """
    get an integer response from the prompt
    """
    if prompt == None:
        prompt = "Enter an integer"
    if default != None:
        if type(default) != types.IntType:
            default = int(default)
        prompt = "%s (%i)" % (prompt,default)
    prompt = prompt + ">>"
    if nlprompt==True: prompt = prompt + "\n"
    while 1:
        try:
            ok = True
            ret = raw_input(prompt)
            if len(ret.strip()) == 0 and default != None:
                return default
            ret = int(ret)
            if len(valid) > 0 and ok == True:
                if ret not in valid:
                    print "Enter integer in range %s " % str(valid)
                    ok = False
            if min != None and ok == True:
                if ret < min:
                    print "Enter integer greater than %i " % min
                    ok = False
            if max != None and ok == True:
                if ret > max:
                    print "Enter integer less than %i " % max
                    ok = False
            if ok == True:
                return ret
        except:
            pass

##########################################################################
def get_flt(prompt=None,default=None,min=None,max=None,nlprompt=True):
    """
    get a float response from the prompt
    """
    if prompt == None:
        prompt = "Enter a float"
    if default != None:
        if type(default) != types.FloatType:
            default = float(default)
        prompt = "%s (%g)" % (prompt,default)
    prompt = prompt + ">>"
    if nlprompt==True: prompt = prompt + "\n"
    while 1:
        try:
            ok = True
            ret = raw_input(prompt)
            if len(ret.strip()) == 0 and default != None:
                return default
            ret = float(ret)
            if min != None and ok == True:
                if ret < min:
                    print "Enter float greater than %g " % min
                    ok = False
            if max != None and ok == True:
                if ret > max:
                    print "Enter float less than %g " % max
                    ok = False
            if ok == True:
                return ret
        except:
            pass

##########################################################################
def get_flt_list(prompt=None,default=None,nlprompt=True):
    """
    get a list of floats from the prompt
    """
    if prompt == None:
        prompt = "Enter a list of floats"
    if default != None:
        prompt = "%s (%s)" % (prompt,str(default))
    prompt = prompt + ">>"
    if nlprompt==True: prompt = prompt + "\n"
    while 1:
        try:
            ret = raw_input(prompt)
            ret = ret.strip()
            if len(ret) == 0:
                if default != None:
                    return default
            else:
                ret = str2list(ret)
                return ret
        except:
            pass

##########################################################################
class Menu:
    """
    A menu driven interface

    Given a set of menu labels
    and corresponding set of descriptions,
    loop until there is a match to one of the labels
    at the command line
    """
    ##########################################################
    def __init__(self,labels=[],descr=[],header=None,
                 dohelp=True,matchidx=True,sort=True):
        """init"""
        self.dohelp = dohelp
        if (len(descr) != 0) and (len(descr) != len(labels)):
            print "Error length mis-match between labels and descriptions"
            raise exceptions.IndexError
        else:
            self._init(labels,descr,header,matchidx,sort)

    ##########################################################
    def _init(self,labels,descr,header,matchidx,sort):
        """init"""
        import numpy as num

        self.matchidx = matchidx

        # add help
        if self.dohelp:
            if 'help' not in labels:
                labels.append('help')
                if len(descr) > 0:
                    descr.append('Show menu options')

        # add labels and descr
        self.labels  = []
        self.descr  = []
        for j in range(len(labels)):
            item = labels[j].lower
            if item in self.labels:
                print "Error '%s' is a repeated menu item" % item
                raise exceptions.IndexError
            else:
                self.labels.append(labels[j])
                if len(descr) > 0:
                    self.descr.append(descr[j])
        nlabels = len(self.labels)

        # sort
        if sort:
            self.labels = num.array(self.labels)
            idx        = num.argsort(self.labels)
            self.labels = list(self.labels[idx][:])
            if len(descr) > 0:
                self.descr = num.array(self.descr)
                self.descr = list(self.descr[idx][:])

        # calc unique array
        # loop through chars of each item and determine
        # when it becomes unique relative to all others
        self.unique  = []
        maxchar = 1
        for j in range(nlabels):
            unique = 1
            nchar  = len(self.labels[j])
            if nchar > maxchar: maxchar = nchar
            for k in range(nchar):
                for l in range(nlabels):
                    if l != j:
                        item = self.labels[j][0:k+1]
                        test = self.labels[l]
                        if test[0:k+1] == item: unique = k+2
            self.unique.append(unique)

        # header
        if header != None:
            self.header = str(header)
        else:
            self.header = None

        # calc options string
        options = ''
        if self.matchidx:
            fmt = "%%s(%%2i) %%-%is: %%s\n" % maxchar
        else:
            fmt = "%%s%%-%is: %%s\n" % maxchar
        for j in range(nlabels):
            if len(self.descr) > 0:
                descr = self.descr[j]
            else:
                descr = ''
            idx  = self.unique[j]
            item = self.labels[j][0:idx].upper() + self.labels[j][idx:]
            if self.matchidx:
                options = fmt % (options,j+1, item, descr)
            else:
                options = fmt % (options, item, descr)
        self.options = options

    ##########################################################
    def match(self,cmd):
        """
        match cmd to list of labels
        """
        cmd0 = cmd
        cmd = cmd.strip().lower()
        nmatch = -1
        cmdnum = -1
        nlabels = len(self.labels)
        cmdlen = len(cmd)

        if self.matchidx:
            try:
                cmdnum = int(cmd)
                if cmdnum in range(1,nlabels+1):
                    cmdnum = cmdnum-1
                    nmatch = 0
            except:
                pass
        if cmdnum == -1:
            for j in range(nlabels):
                if (cmdlen >= self.unique[j]):
                    if self.labels[j].find(cmd,0,cmdlen) == 0:
                        nmatch = nmatch + 1
                        cmdnum = j
        if nmatch == -1:
            print "Command '%s' not found" % cmd0
            return None
        elif nmatch > 0:
            print "Warning: '%s' is an ambiguous command" % cmd0
            return None
        else:
            if self.dohelp == True:
                if self.labels[cmdnum] == 'help':
                    print self.options
                    return None
            return self.labels[cmdnum]

    ##########################################################
    def prompt(self,p='>>',nlprompt=True):
        """
        Issue prompt till we recieve a match
        Should we add a default option?
        """
        if self.header != None: print self.header
        if self.options != None: print self.options
        #
        def _input(p):
            cmd = raw_input(p)
            cmd = cmd.strip()
            if len(cmd) == 0:
                return None
            else:
                return cmd
        #
        if nlprompt==True:
            p = p + "\n"
        cmd = _input(p)
        if cmd != None:
            cmd = self.match(cmd)
        while cmd == None:
            cmd = _input(p)
            if cmd != None:
                cmd = self.match(cmd)
        return cmd

##########################################################################
##########################################################################
if __name__ == "__main__":
    labels = ['save','quit','run', 'zzz','ran','runner','n']
    descr = ['save stuff',
             'all done',
             'do a thing',
             'zzzzzzz',
             'do the other thing',
             'x','y']
    header = "adf\njkd\n"
    m = Menu(labels=labels,descr=descr,header=header,sort=False,matchidx=False)
    print m.prompt('xx>')
    #
    print get_yn()
    print get_tf()
    print get_int()

