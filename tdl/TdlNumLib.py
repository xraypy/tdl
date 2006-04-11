# M. Newville Univ of Chicago (2005)
#
# --------------
# Modifications
# --------------
#
# * 1-28-06 T2:
#  - Created, moved numeric related funs/cmds to this module
#   from TdlLib
##########################################################################
from Num import Num
import os
import sys
import types


title = "numeric functions"

def _identity(x):
    "wrap Numeric.identity"
    return Num.identity(int(x))

def _take(x,y):
    "wrap Numeric.take"
    return Num.take(x,y.astype(Num.Int32))

def _choose(x,y):
    "wrap Numeric.choose"
    return Num.choose(x.astype(Num.Int32),y)

def _int(x):
    "wrap builtin int or Numpy astype(int)"
    if type(x) == Num.ArrayType:
        return x.astype(int)
    else:
        return int(x)

def _float(x):
    "wrap builtin float or numpy astype(float)"
    if type(x) == Num.ArrayType:
        return x.astype(float)
    else:
        return float(x)

def _complex(x):
    "wrap builtin complex or numpy astype(complex)"
    if type(x) == Num.ArrayType:
        return x.astype(complex)
    else:
        return complex(x)


def _range(x,stop=None,step=None,shape=None,dtype='d'):
    """create an array of evenly spaced values:
        range(x, [stop=stop, [step=step, [shape=shape, [dtype='d']]]])

    Thus,
       range(10)
    returns [0. 1. 2. 3. 4. 5. 6. 7. 8. 9.] (the first 10 numers, as double precision)
    other variations: 
       range(2,5)               -> [2. 3. 4.]
       range(2,10,2)            -> [2. 4. 6. 8.]
       range(2,10,2,dtype='i')  -> [2 4 6 8]

       range(20,shape=[5,4])    ->  [[  0.   1.   2.   3.]
                                     [  4.   5.   6.   7.]
                                     [  8.   9.  10.  11.]
                                     [ 12.  13.  14.  15.]
                                     [ 16.  17.  18.  19.]]
    """
    if stop is None and step is None:             t = Num.arange(x,dtype=dtype)
    elif stop is not None and step is None:       t = Num.arange(x,stop,dtype=dtype)
    elif stop is not None and step is not None:   t = Num.arange(x,stop,step,dtype=dtype)
    elif stop is None     and step is not None:   t = Num.arange(0, x,step,dtype=dtype)

    if shape is None: return t
    reshaped = False
    if type(shape) in  (types.ListType,types.TupleType):
        npts = 1
        for i in shape: npts = npts*int(i)
        if npts == len(t):
            try:
                t = Num.array(t)
                t.shape = tuple(shape)
                return t
            except ValueError:
                pass
    # if we got here, there's a value error
    raise ValueError, ' could not re-shape array to shape = %s ' % repr(shape)
            

#################################################################
# Load the functions
#################################################################

# constants to go into _builtin name space
_consts_ = {"_math.pi":Num.pi, "_math.e":Num.e}


# these do not get passed a reference to tdl in call argument
# should we make all these num.fcn ???
_func_ = {"_math.array":Num.array,
          "_math.sin":Num.sin,
          "_math.cos":Num.cos,
          "_math.tan":Num.tan,          
          "_math.exp":Num.exp,
          "_math.ln":Num.log,
          "_math.log":Num.log,
          "_math.log10":Num.log10,
          "_math.sqrt":Num.sqrt,
          "_math.int": _int,
          "_math.float": _float,
          "_math.complex": _complex,
          "_math.range":_range,
          "_math.arange":_range,
          "_math.arccos":Num.arccos,
          "_math.arccosh":Num.arccosh,
          "_math.arcsin":Num.arcsin,
          "_math.arcsinh":Num.arcsinh,
          "_math.arctan":Num.arctan,
          "_math.arctan2":Num.arctan2,
          "_math.arctanh":Num.arctanh,
          "_math.argmax":Num.argmax,
          "_math.argmin":Num.argmin,
          "_math.argsort":Num.argsort,
          "_math.array":Num.array,
          "_math.arrayrange":Num.arrayrange,
          "_math.cosh":Num.cosh,
          "_math.exp":Num.exp,
          "_math.fabs":Num.fabs,
          "_math.floor":Num.floor,
          "_math.floor_divide":Num.floor_divide,
          "_math.fmod":Num.fmod,
          "_math.log":Num.log,
          "_math.log10":Num.log10,
          "_math.tan":Num.tan,
          "_math.tanh":Num.tanh,
          "_math.sign":Num.sign,
          "_math.sin":Num.sin,
          "_math.sinh":Num.sinh,
          "_math.sqrt":Num.sqrt,          
          "array.identity": _identity,
          "array.take": _take,
          "array.choose": _choose,
          "array.add":Num.add,
          "array.allclose":Num.allclose,
          "array.alltrue":Num.alltrue,
          "array.around":Num.around,
          "array.asarray":Num.asarray,
          "array.average":Num.average,
          "array.bitwise_and":Num.bitwise_and,
          "array.bitwise_or":Num.bitwise_or,
          "array.bitwise_xor":Num.bitwise_xor,
          "array.ceil":Num.ceil,
          "array.clip":Num.clip,
          "array.compress":Num.compress,
          "array.concatenate":Num.concatenate,
          "array.conjugate":Num.conjugate,
          "array.convolve":Num.convolve,
          "array.cross_correlate":Num.cross_correlate,
          "array.cumproduct":Num.cumproduct,
          "array.cumsum":Num.cumsum,
          "array.diagonal":Num.diagonal,
          "array.divide":Num.divide,
          "array.dot":Num.dot,
          "array.dump":Num.dump,
          "array.dumps":Num.dumps,
          "array.equal":Num.equal,
          "array.fromfunction":Num.fromfunction,
          "array.fromstring":Num.fromstring,
          "array.greater":Num.greater,
          "array.greater_equal":Num.greater_equal,
          "array.hypot":Num.hypot,
          "array.indices":Num.indices,
          "array.innerproduct":Num.innerproduct,
          "array.invert":Num.invert,
          "array.left_shift":Num.left_shift,
          "array.less":Num.less,
          "array.less_equal":Num.less_equal,
          "array.load":Num.load,
          "array.loads":Num.loads,
          "array.logical_and":Num.logical_and,
          "array.logical_not":Num.logical_not,
          "array.logical_or":Num.logical_or,
          "array.logical_xor":Num.logical_xor,
          "array.matrixmultiply":Num.matrixmultiply,
          "array.maximum":Num.maximum,
          "array.minimum":Num.minimum,
          "array.multiply":Num.multiply,
          "array.negative":Num.negative,
          "array.nonzero":Num.nonzero,
          "array.not_equal":Num.not_equal,
          "array.ones":Num.ones,
          "array.outerproduct":Num.outerproduct,
          "array.power":Num.power,
          "array.product":Num.product,
          "array.put":Num.put,
          "array.putmask":Num.putmask,
          "array.rank":Num.rank,
          "array.ravel":Num.ravel,
          "array.remainder":Num.remainder,
          "array.repeat":Num.repeat,
          "array.reshape":Num.reshape,
          "array.resize":Num.resize,
          "array.right_shift":Num.right_shift,
          "array.searchsorted":Num.searchsorted,
          "array.shape":Num.shape,
          "array.size":Num.size,
          "array.sometrue":Num.sometrue,
          "array.sort":Num.sort,
          "array.subtract":Num.subtract,
          "array.sum":Num.sum,
          "array.swapaxes":Num.swapaxes,
          "array.trace":Num.trace,
          "array.transpose":Num.transpose,
          "array.true_divide":Num.true_divide,
          "array.vdot":Num.vdot,
          "array.where":Num.where,
          "array.zeros":Num.zeros,
          }
          


if __name__ == '__main__':
    print help
