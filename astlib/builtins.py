import __builtin__
import copy
import numpy
from util import closure

# inherit these from python's __builtin__
_from_builtin= ('ArithmeticError', 'AssertionError', 'AttributeError',
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
               'true_divide', 'vdot', 'where', 'zeros','linspace')

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
        
def load_all(compiler):
    symtable      = compiler.symtable
    import_module = compiler.import_module

    imports = ((_from_builtin, __builtin__,'_builtin'),
               (_from_numpy,   numpy, '_math'))
               
    for symlist, pymod, tdlmod in imports:

        group = getattr(symtable,tdlmod)
        for sym in symlist:
            setattr(group,sym,getattr(pymod,sym))

    group = getattr(symtable,'_builtin')
    for fname,fcn in _local_funcs.items():
        setattr(group, fname, closure(func=fcn))
