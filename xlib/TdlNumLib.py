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
from Num import Num, num_version
import os
import sys
import types
from Util import datalen, EvalException

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

def _min(x,*args):
    "return mininum value of an array or list"
    t = []
    for i in args + (x,):
        if type(i) == Num.ArrayType: i = i.min()
        t.append(i)
    return min(t)


def _max(x,*args):
    "return maxinum value of an array or list"    
    t = []
    for i in args + (x,):
        if type(i) == Num.ArrayType: i = i.max()
        t.append(i)
    return max(t)

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
            

def _tdl_minimize(residual,varlist,tdl=None, iter_function=None, debug=False, **kw):
    """tdl optimization:
         minimize(f, variable_list)
    """
    if tdl is None: return varlist
    if num_version.startswith('numpy'):
        print 'cannot use tdl minimize:  scipy required'
        return varlist

    from scipy.optimize import leastsq
    tdl.setVariable('_sys.fit_iterations',0)

    params = []
    for i in varlist: params.append(tdl.getVariableValue(i))
    params = Num.array(params)

    tdl.eval(residual)

    rsym = tdl.getVariable(residual)
    try:
        if rsym.type in ('pyfunc','defvar'): pass
    except:
        print 'cannot minimize %s: must be a defined procedure or variable ' % residual
        return varlist
    
    n_iter = tdl.getVariable('_sys.fit_iterations')
    def get_residual(params):
        for name,val in zip(varlist,params): tdl.setVariable(name,val)
        r = tdl.eval(residual)
        n_iter.value = n_iter.value + 1
        if iter_function is not None:  tdl.eval(iter_function)
        return r

    result = leastsq(get_residual,params)
    tdl.eval(residual)
    tdl.setVariable('_sys.fit_message',result[1])
    if debug: print  ' optimization done:  ', n_iter.value, result[1]
    return result[0]

def _random_seed(x=None):
    "wrap numpy random seed "
    ## print x, datalen(x)
    if x is None:
        return Num.random.seed()
    else:
        try:
            return Num.random.seed([x])
        except:
            return Num.random.seed()            

def _random(a,b=1,c=1,npts=1,distribution='normal',**kw):
    "wrapper for numpy random distributions" 
    NR = Num.random
    if   distribution == 'binomial':        return NR.binomial(a,b,size=npts)
    elif distribution == 'geometric':       return NR.geometric(a,size=npts)    
    elif distribution == 'poisson':         return NR.poisson(a,size=npts)    
    elif distribution == 'zipf':            return NR.zipf(a,size=npts)    
    elif distribution == 'beta':            return NR.beta(a,b,size=npts)    
    elif distribution == 'chisquare':       return NR.chisquare(a,size=npts)    
    elif distribution == 'exponential':     return NR.exponential(a,size=npts)
    elif distribution == 'gamma':           return NR.gamma(a,b,size=npts)
    elif distribution == 'gumbel':          return NR.gumbel(a,b,size=npts)
    elif distribution == 'laplace':         return NR.laplace(a,b,size=npts)
    elif distribution == 'lognormal':       return NR.lognormal(a,b,size=npts)    
    elif distribution == 'logistic':        return NR.logistic(a,b,size=npts)    
    elif distribution == 'multivariate_normal':   return NR.multivariate_normal(a,b,size=npts)
    elif distribution == 'noncentral_chisquare':  return NR.noncentral_chisquare(a,b,size=npts)
    elif distribution == 'noncentral_f':    return NR.noncentral_f(a,b,c,size=npts)
    elif distribution == 'normal':          return NR.normal(a,b,size=npts)
    elif distribution == 'pareto':          return NR.pareto(a,size=npts)
    elif distribution == 'power':           return NR.power(a,size=npts)
    elif distribution == 'randint':         return NR.randint(a,b,size=npts)
    elif distribution == 'random_integers': return NR.random_integers(a,b,size=npts)
    elif distribution == 'rayleigh':        return NR.rayleigh(a,size=npts)
    elif distribution == 'standard_cauchy': return NR.standard_cauchy(size=npts)
    elif distribution == 'standard_exponential':        return NR.standard_exponential(size=npts)
    elif distribution == 'standard_gamma':  return NR.standard_gamma(a,size=npts)
    elif distribution == 'standard_normal': return NR.standard_normal(size=npts)
    elif distribution == 'standard_t':      return NR.standard_t(a,size=npts)
    elif distribution == 'uniform':         return NR.uniform(a,b,size=npts)
    elif distribution == 'wald':            return NR.wald(a,b,size=npts)
    elif distribution == 'weibull':         return NR.weibull(a,b,size=npts)


#################################################################
# Load the functions
#################################################################

# constants to go into _builtin name space
_consts_ = {"_math.pi":Num.pi, "_math.e":Num.e}


# these do not get passed a reference to tdl in call argument
# should we make all these num.fcn ???
_func_ = {
          "_math.minimize":_tdl_minimize,
          "_math.array":Num.array,
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
          "_math.random_seed": _random_seed,
          "_math.random": _random,
          "_math.identity": _identity,
          "_math.take": _take,
          "_math.choose": _choose,
          "_math.add":Num.add,
          "_math.allclose":Num.allclose,
          "_math.alltrue":Num.alltrue,
          "_math.around":Num.around,
          "_math.asarray":Num.asarray,
          "_math.average":Num.average,
          "_math.bitwise_and":Num.bitwise_and,
          "_math.bitwise_or":Num.bitwise_or,
          "_math.bitwise_xor":Num.bitwise_xor,
          "_math.ceil":Num.ceil,
          "_math.clip":Num.clip,
          "_math.compress":Num.compress,
          "_math.concatenate":Num.concatenate,
          "_math.conjugate":Num.conjugate,
          "_math.convolve":Num.convolve,
          "_math.cross_correlate":Num.cross_correlate,
          "_math.cumproduct":Num.cumproduct,
          "_math.cumsum":Num.cumsum,
          "_math.diagonal":Num.diagonal,
          "_math.divide":Num.divide,
          "_math.dot":Num.dot,
          "_math.dump":Num.dump,
          "_math.dumps":Num.dumps,
          "_math.equal":Num.equal,
          "_math.fromfunction":Num.fromfunction,
          "_math.fromstring":Num.fromstring,
          "_math.greater":Num.greater,
          "_math.greater_equal":Num.greater_equal,
          "_math.hypot":Num.hypot,
          "_math.indices":Num.indices,
          "_math.innerproduct":Num.innerproduct,
          "_math.invert":Num.invert,
          "_math.left_shift":Num.left_shift,
          "_math.less":Num.less,
          "_math.less_equal":Num.less_equal,
          # "_math.load":Num.load,
          # "_math.loads":Num.loads,
          "_math.logical_and":Num.logical_and,
          "_math.logical_not":Num.logical_not,
          "_math.logical_or":Num.logical_or,
          "_math.logical_xor":Num.logical_xor,
          "_math.matrixinverse":Num.linalg.inv,
          "_math.matrixmultiply":Num.matrixmultiply,
          "_math.maximum":Num.maximum,
          "_math.minimum":Num.minimum,
          "_math.multiply":Num.multiply,
          "_math.negative":Num.negative,
          "_math.nonzero":Num.nonzero,
          "_math.not_equal":Num.not_equal,
          "_math.ones":Num.ones,
          "_math.outerproduct":Num.outerproduct,
          "_math.power":Num.power,
          "_math.product":Num.product,
          "_math.put":Num.put,
          "_math.putmask":Num.putmask,
          "_math.rank":Num.rank,
          "_math.ravel":Num.ravel,
          "_math.remainder":Num.remainder,
          "_math.repeat":Num.repeat,
          "_math.reshape":Num.reshape,
          "_math.resize":Num.resize,
          "_math.right_shift":Num.right_shift,
          "_math.searchsorted":Num.searchsorted,
          "_math.shape":Num.shape,
          "_math.size":Num.size,
          "_math.sometrue":Num.sometrue,
          "_math.sort":Num.sort,
          "_math.subtract":Num.subtract,
          "_math.sum":Num.sum,
          "_math.swapaxes":Num.swapaxes,
          "_math.trace":Num.trace,
          "_math.transpose":Num.transpose,
          "_math.true_divide":Num.true_divide,
          "_math.vdot":Num.vdot,
          "_math.where":Num.where,
          "_math.zeros":Num.zeros,
          "_builtin.max":_max,
          "_builtin.min":_min,
          }
          

if __name__ == '__main__':
    print help
