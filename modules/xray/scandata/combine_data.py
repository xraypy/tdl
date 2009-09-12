#######################################################################
"""
T. Trainor (fftpt@uaf.edu)
Functions for averaging and merging ScanData objects

Modifications:
--------------

"""
#######################################################################
"""
Todo

- Test!

"""
########################################################################

import types
import copy
import numpy as num

from scandata import ScanData

########################################################################
def append_data(data1,data2,sort=True):
    """
    Append 2 ScanData instances
    This works for one dimensional data

    Note check how this works.  only append arrays of full length
    others combine into 2d arrays.  is this ok???  what if combine three
    or more... 
    
    """
    # Check data
    if (len(data1.dims) != 1) or (len(data1.primary_axis) > 1):
        print "Only 1-d data sets may be combined"
        return
    if len(data2.dims) != 1 or (len(data2.primary_axis) > 1):
        print "Only 1-d data sets may be combined"
        return
    if data2.primary_axis != data1.primary_axis:
        print "Warning primary axis doesnt match "
    if data2.primary_det != data1.primary_det:
        print "Warning primary detector doesnt match "
    if data1.xrf_lines != data2.xrf_lines:
        print "Warning xrf lines dont match"
    if data1.image_rois != data2.image_rois:
        print "Warning image rois dont match"

    # init
    data3 = ScanData()
    data3.name         = "%s, %s" % (data1.name, data2.name)
    data3.dims         = [data1.dims[0]+data2.dims[0]]
    npts               = data3.dims[0]
    data3.primary_axis = copy.copy(data1.primary_axis)
    data3.primary_det  = copy.copy(data1.primary_det)

    # append scalers, all will be numarrays
    for key in data1.scalers.keys():
        s1 = num.array(data1.scalers[key])
        s2 = num.array(data2.scalers[key])
        #tmp = num.concatenate( (s1,s2) )
        tmp = num.append( s1,s2 )
        if len(tmp) == npts:
            data3.scalers.update({key:tmp})
        else:
            data3.scalers.update({key:num.array([s1,s2])})

    # append positioners, all will be numarrays
    for key in data1.positioners.keys():
        s1 = num.array(data1.positioners[key])
        s2 = num.array(data2.positioners[key])
        tmp = num.append(s1,s2)
        if len(tmp) == npts:
            data3.positioners.update({key:tmp})
        else:
            data3.positioners.update({key:num.array([s1,s2])})

    # combine state info,
    # these are not appened!  
    for key in data1.state.keys():
        s1 = num.array(data1.state[key])
        s2 = num.array(data2.state[key])
        tmp = num.array([s1,s2])
        data3.state.update({key:tmp})

    # append meds, xrfs...
    data3.med = num.append(data1.med, data2.med)
    data3.xrf = num.append(data1.xrf, data2.xrf)
    data3.xrf_lines = data1.xrf_lines
    for key in data1.xrf_peaks.keys():
        s1 = num.array(data1.xrf_peaks[key])
        s2 = num.array(data2.xrf_peaks[key])
        tmp = num.append(s1,s2)
        data3.xrf_peaks.update({key:tmp})
    
    # append images...
    data3.image = data1.image
    for i in range(data2.dims[0]):
        data3.image.append(data2.image[i])

    data3.image_rois = data1.image_rois
    for key in data1.image_peaks.keys():
        s1 = num.array(data1.image_peaks[key])
        s2 = num.array(data2.image_peaks[key])
        tmp = num.append(s1,s2)
        data3.image_peaks.update({key:tmp})

    if sort: _sort_data(data3)
        
    return data3

def _sort_data(data):
    """
    sort 1d data on the principle axis
    """
    paxis = data.primary_axis
    if len(paxis) != 1: return
    
    a    = data.positioners[paxis[0]]
    npts = len(data.positioners[paxis[0]])
    idx  = a.argsort()

    for key in data.scalers.keys():
        if len(data.scalers[key]) == npts:
            data.scalers[key] = data.scalers[key][idx]
    
    for key in data.positioners.keys():
        if len(data.positioners[key]) == npts:
            data.positioners[key] = data.positioners[key][idx]

    if len(data.med) == npts:
        data.med = data.med[idx]

    if len(data.xrf) == npts:
        data.xrf = data.xrf[idx]

    for key in data.xrf_peaks.keys():
        if len(data.xrf_peaks[key]) == npts:
            data.xrf_peaks[key] = data.xrf_peaks[key][idx]

    if len(data.image) == npts:
        data.image = data.image[idx]

    for key in data.image_peaks.keys():
        if len(data.image_peaks[key]) == npts:
            data.image_peaks[key] = data.image_peaks[key][idx]

    return            

########################################################################
def merge_data(data=[],average=True,align=False,fast=True):
    """
    Merge  a list of ScanData instances
    This works for one dimensional data
    
    """
    ndat = len(data)
    if ndat < 2: return None
    
    # The first one in the list sets the dimensions and
    # default values etc... check that all agree
    name = ''
    for d in data:
        if (len(d.dims) != 1) or (len(d.primary_axis) > 1):
            print "Only 1-d data sets may be combined"
            return
        if d.primary_axis != data[0].primary_axis:
            print "Warning primary axis doesnt match "
        if d.primary_det != data[0].primary_det:
            print "Warning primary detector doesnt match "
        if d.xrf_lines != data[0].xrf_lines:
            print "Warning xrf lines dont match"
        if d.image_rois != data[0].image_rois:
            print "Warning image rois dont match"
        if len(name) == 0:
            name = "%s" % d.name
        else:
            name = "%s, %s" % (name, d.name)
        
    # init
    data_m      = ScanData()
    data_m.name = name
    data_m.dims = data[0].dims
    npts        = data[0].dims[0]
    paxis = data_m.primary_axis = copy.copy(data[0].primary_axis)
    pdet  = data_m.primary_det  = copy.copy(data[0].primary_det)

    # the primary axis of data[0] sets the primary
    # axis of data_m when align is true
    if align:
        newx = data[0][paxis]
    else:
        newx = None

    # sum/average scalers
    for key in data[0].scalers.keys():
        tmp = []
        for j in range(ndat):
            if j == 0:
                tmp.append(num.array(data[j].scalers[key]))
            else:
                s = num.array(data[j].scalers[key])
                if align and (len(s) == npts):
                    oldx = data[j][paxis]
                    s = spline_interpolate(oldx,s,newx,fast=fast)
                tmp.append(s)
        if len(tmp[0]) == npts:
            tmp = num.sum(tmp,0)
            if average:
                tmp = tmp / (1.0*ndat)
        data_m.scalers.update({key:tmp})


    # average positioners
    for key in data[0].positioners.keys():
        tmp = []
        for j in range(ndat):
            if j == 0:
                tmp.append(num.array(data[j].positioners[key]))
            else:
                s = num.array(data[j].positioner[key])
                if align and (len(s) == npts):
                    # does this work for positioners?
                    oldx = data[j][paxis]
                    s = spline_interpolate(oldx,s,newx,fast=fast)
                tmp.append(s)
        if len(tmp[0]) == npts:
            tmp = num.sum(tmp,0)
            # always average
            tmp = tmp / (1.0*ndat)
        data_m.scalers.update({key:tmp})

    # combine state info,
    # these are not summed/averaged!  
    for key in data[0].state.keys():
        tmp = []
        for j in range(ndat):
            tmp.append(data[j].state[key])
        data_m.state.update({key:tmp})

    # append meds, xrfs, images...
    data_m.xrf_lines = data[0].xrf_lines
    data_m.image_rois = data[0].image_rois

    # should we try to sum the spectra (and images)?
    for j in range(ndat):
        data_m.med.append(data[j].med)
        data_m.xrf.append(data[j].xrf)
        data_m.image.append(data[j].image)
        

    # treat peaks etc same as scalers....
    for key in data[0].xrf_peaks.keys():
        tmp = []
        for j in range(ndat):
            if j == 0:
                tmp.append(num.array(data[j].xrf_peaks[key]))
            else:
                s = num.array(data[j].xrf_peaks[key])
                if align and (len(s) == npts):
                    oldx = data[j][paxis]
                    s = spline_interpolate(oldx,s,newx,fast=fast)
                tmp.append(s)
        if len(tmp[0]) == npts:
            tmp = num.sum(tmp,0)
            if average:
                tmp = tmp / (1.0*ndat)
        data_m.xrf_peaks.update({key:tmp})

    for key in data[0].image_peaks.keys():
        tmp = []
        for j in range(ndat):
            if j == 0:
                tmp.append(num.array(data[j].image_peaks[key]))
            else:
                s = num.array(data[j].image_peaks[key])
                if align and (len(s) == npts):
                    oldx = data[j][paxis]
                    s = spline_interpolate(oldx,s,newx,fast=fast)
                tmp.append(s)
        if len(tmp[0]) == npts:
            tmp = num.sum(tmp,0)
            if average:
                tmp = tmp / (1.0*ndat)
        data_m.image_peaks.update({key:tmp})
        
    return data_m

################################################################################
################################################################################
"""
cubic splines for axis alignment using
scipy.signal and/or scipy.interpolate
"""
from scipy.interpolate import splrep, splev
from scipy.signal import cspline1d, cspline1d_eval

def spline_interpolate(oldx, oldy, newx, smoothing=0.001,fast=True, **kw):
    """
    newy = spline_interpolate(oldx, oldy, newx, fast=True)
    if fast = True
       1-dimensional cubic spline for cases where
       oldx and newx are on a uniform grid.
    else
       handles multi-dimensional data, non-uniform x-grids, but is
       much slower for 1d cubic splines
    """
    if fast:
        return cspline1d_eval(cspline1d(oldy), newx, dx=oldx[1]-oldx[0],x0=oldx[0])
    else:
        rep = splrep(oldx,oldy,s=smoothing,full_output=False,**kw)
        return splev(newx, rep)

################################################################################

