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

#################################################################
# Load the functions
#################################################################

# constants to go into _builtin name space
# _consts_ = {"pi":Num.pi, "e":Num.e}

# these do not get passed a reference to tdl in call argument
# should we make all these num.fcn ???
_func_ = {"num.identity": _identity,
          "num.take": _take,
          "num.choose": _choose,
          "num.arange":Num.arange,
          "num.add":Num.add,
          "num.allclose":Num.allclose,
          "num.alltrue":Num.alltrue,
          "num.arccos":Num.arccos,
          "num.arccosh":Num.arccosh,
          "num.arcsin":Num.arcsin,
          "num.arcsinh":Num.arcsinh,
          "num.arctan":Num.arctan,
          "num.arctan2":Num.arctan2,
          "num.arctanh":Num.arctanh,
          "num.argmax":Num.argmax,
          "num.argmin":Num.argmin,
          "num.argsort":Num.argsort,
          "num.around":Num.around,
          "num.array":Num.array,
          "num.arrayrange":Num.arrayrange,
          "num.asarray":Num.asarray,
          "num.average":Num.average,
          "num.bitwise_and":Num.bitwise_and,
          "num.bitwise_or":Num.bitwise_or,
          "num.bitwise_xor":Num.bitwise_xor,
          "num.ceil":Num.ceil,
          "num.clip":Num.clip,
          "num.compress":Num.compress,
          "num.concatenate":Num.concatenate,
          "num.conjugate":Num.conjugate,
          "num.convolve":Num.convolve,
          "num.cos":Num.cos,
          "num.cosh":Num.cosh,
          "num.cross_correlate":Num.cross_correlate,
          "num.cumproduct":Num.cumproduct,
          "num.cumsum":Num.cumsum,
          "num.diagonal":Num.diagonal,
          "num.divide":Num.divide,
          "num.dot":Num.dot,
          "num.dump":Num.dump,
          "num.dumps":Num.dumps,
          "num.equal":Num.equal,
          "num.exp":Num.exp,
          "num.fabs":Num.fabs,
          "num.floor":Num.floor,
          "num.floor_divide":Num.floor_divide,
          "num.fmod":Num.fmod,
          "num.fromfunction":Num.fromfunction,
          "num.fromstring":Num.fromstring,
          "num.greater":Num.greater,
          "num.greater_equal":Num.greater_equal,
          "num.hypot":Num.hypot,
          "num.indices":Num.indices,
          "num.innerproduct":Num.innerproduct,
          "num.invert":Num.invert,
          "num.left_shift":Num.left_shift,
          "num.less":Num.less,
          "num.less_equal":Num.less_equal,
          "num.load":Num.load,
          "num.loads":Num.loads,
          "num.log":Num.log,
          "num.log10":Num.log10,
          "num.logical_and":Num.logical_and,
          "num.logical_not":Num.logical_not,
          "num.logical_or":Num.logical_or,
          "num.logical_xor":Num.logical_xor,
          "num.matrixmultiply":Num.matrixmultiply,
          "num.maximum":Num.maximum,
          "num.minimum":Num.minimum,
          "num.multiply":Num.multiply,
          "num.negative":Num.negative,
          "num.nonzero":Num.nonzero,
          "num.not_equal":Num.not_equal,
          "num.ones":Num.ones,
          "num.outerproduct":Num.outerproduct,
          "num.power":Num.power,
          "num.product":Num.product,
          "num.put":Num.put,
          "num.putmask":Num.putmask,
          "num.rank":Num.rank,
          "num.ravel":Num.ravel,
          "num.remainder":Num.remainder,
          "num.repeat":Num.repeat,
          "num.reshape":Num.reshape,
          "num.resize":Num.resize,
          "num.right_shift":Num.right_shift,
          "num.searchsorted":Num.searchsorted,
          "num.shape":Num.shape,
          "num.sign":Num.sign,
          "num.sin":Num.sin,
          "num.sinh":Num.sinh,
          "num.size":Num.size,
          "num.sometrue":Num.sometrue,
          "num.sort":Num.sort,
          "num.sqrt":Num.sqrt,
          "num.subtract":Num.subtract,
          "num.sum":Num.sum,
          "num.swapaxes":Num.swapaxes,
          "num.tan":Num.tan,
          "num.tanh":Num.tanh,
          "num.trace":Num.trace,
          "num.transpose":Num.transpose,
          "num.true_divide":Num.true_divide,
          "num.vdot":Num.vdot,
          "num.where":Num.where,
          "num.zeros":Num.zeros}


if __name__ == '__main__':
    print help
