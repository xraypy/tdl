""" cubic spline:  now using scipy.interpolate splrep and splev """

from scipy.interpolate import splrep, splev
def spline_interpolate(oldx, oldy, newx, smoothing=0.001, **kw):
    """ newy = spline_interpolate(oldx, oldy, newx)

    """
    rep = splrep(oldx,oldy,s=smoothness,full_output=False,**kw)
    return splev(newx, rep)
