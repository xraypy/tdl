""" cubic spline:  now using scipy.signal and/or scipy.interpolate """

from scipy.interpolate import splrep, splev
from scipy.signal import cspline1d, cspline1d_eval

def spline_interpolate_general(oldx, oldy, newx, smoothing=0.001, **kw):
    """ newy = spline_interpolate_general(oldx, oldy, newx)

     handles multi-dimensional data, non-uniform x-grids, but is
     much slower than spline_interpolate for 1d cubic splines

    """
    rep = splrep(oldx,oldy,s=smoothing,full_output=False,**kw)
    return splev(newx, rep)

def spline_interpolate(oldx, oldy, newx, smoothing=0.001, **kw):
    """
      newy = spline_interpolate(oldx, oldy, newx)

      1-dimensional cubic spline, for cases where oldx and newx are on a uniform grid.
    """
    return cspline1d_eval(cspline1d(oldy), newx, dx=oldx[1]-oldx[0],x0=oldx[0])
