#######################################################################
"""
Tom Trainor (fftpt@uaf.edu)
Class to store and operate on generic 'scan data'

Modifications:
--------------

"""
#######################################################################
"""
Todo
- include a flag to control image integration?
- set individual image background parameters for each scan point
- use the bad_points list in image_integrate
"""
########################################################################

import types
from   glob import glob
import os, copy
import numpy as num

import image_data
import med_data
import xrf_data

########################################################################
class ScanData:
    """
    ScanData Class holds:
    name = scan name
    dims ->[npts1] for a 1-d scan
           [npts1,npts2] for a 2-d scan
            where npts1 is the outer loop
            and npts2 is for the inner loop
    scalers = {'io':[],'i1':[]}  -> each array has dims = dims
    positioners  = {'E':[]}      -> these may be single values,
                                    or same size as scalers
    primary_axis = ['E']         -> if multi-dim, list of outer, inner
    primary_det  = 'i1'          -> use for default plotting
    med = []                     -> array has dims = dims or is [] for No med
    xrf = []                     -> array has dims = dims or is [] for No xrf
    images = []                  -> array has dims = dims or is [] for No images
    state = {}                   -> dictionary of additional state information
    """
    ################################################################
    def __init__(self,name='',dims=[],scalers={},positioners={},
                 primary_axis=[],primary_det=None,state={},
                 med=[],xrf=[],xrf_lines=None,
                 image=[],image_rois=[]):
        
        self.name         = name
        self.dims         = dims
        self.primary_axis = primary_axis
        self.primary_det  = primary_det
        self.scalers      = scalers
        self.positioners  = positioners
        self.state        = state
        #
        if len(med)>0:
            self.med = med_data.MedScan(med)
        else:
            self.med = []
        #
        if len(xrf)>0:
            self.xrf = xrf_data.XrfScan(xrf,xrf_lines)
        else:
            self.xrf = []
        #
        #self.xrf          = xrf
        #self.xrf_lines    = xrf_lines
        #self.xrf_peaks    = {}
        #
        if len(image)>0:
            self.image = image_data.ImageScan(image,image_rois)
        else:
            self.image = []
        #
        self.bad_points   = []

    ########################################################################
    """
    def __setitem__(self,arg,val):
        # note this currenlty craps for sub arrays.
        # ie x[0][0].  In such a case arg is a tuple
        # first is the name, the rest are the subscripts
        #print arg, type(arg)
        #print val, type(val)
        if type(arg) == types.StringType:
            if arg.startswith("_"):
                raise TypeError, "Private attribute"
        #if hasattr(self,arg):
        #    setattr(self,arg,val)
        setattr(self,arg,val)
        return
    """
    ########################################################################
    def __getitem__(self,arg):
        """
        Get items.  Note that order matters:
        spectra, images, scalers, positioners, state
        """
        if type(arg) == types.StringType:
            #if arg.startswith("_"):
            #    raise TypeError, "Private attribute"
            #return getattr(self,arg)
            if arg == 'med':
                return self.med
            if arg == 'xrf':
                return self.xrf
            if arg == 'image':
                return self.image
            #
            for key in self.xrf_peaks.keys():
                if key == arg:
                    return self.xrf_peaks[key]
            #
            for key in self.image_peaks.keys():
                if key == arg:
                    return self.image_peaks[key]
            #
            for key in self.scalers.keys():
                if key == arg:
                    return self.scalers[key]
            for key in self.positioners.keys():
                if key == arg:
                    return self.positioners[key]
            for key in self.state.keys():
                if key == arg:
                    return self.state[key]
        return None

    ########################################################################
    def __repr__(self):
        lout = "* Scan Name = %s\n" % self.name
        lout = lout + "* Scan Dimensions = %s\n" % str(self.dims)

        lout = lout + "* Scalers:\n"
        ct = 0
        for sc in self.scalers.keys():
            lout = lout + "%12s" % sc
            ct = ct + 1
            if ct > 7:
                lout = lout + '\n'
                ct = 0

        lout = lout + "\n* Positioners:\n"
        ct = 0
        for p in self.positioners.keys():
            lout = lout + "%12s" % p
            ct = ct + 1
            if ct > 7:
                lout = lout + '\n'
                ct = 0

        if len(self.state) > 0:
            lout = lout + "\n* Additional State Variables:\n"  
            ct = 0
            for p in self.state.keys():
                lout = lout + "%12s" % p
                ct = ct + 1
                if ct > 7:
                    lout = lout + '\n'
                    ct = 0

        lout = lout + "\n* Primary scan axis = %s\n"  % str(self.primary_axis)
        lout = lout + "* Primary detector  = %s\n"  % self.primary_det

        spectra = False
        if self.med != []:
            spectra = True
            lout = lout + "* Scan includes %i med spectra\n" % len(self.med)
        if self.xrf != []:
            spectra = True
            lout = lout + "* Scan includes %i xrf spectra\n" % len(self.xrf)
            if self.xrf_lines != None:
                lout = lout + "  -> Xrf lines = %s\n" % str(self.xrf_lines)
        if spectra == False:
            lout = lout + "* Scan does not include spectra\n"

        if self.image != []:
            lout = lout + "* Scan includes %i images\n" % len(self.image)
        else:
            lout = lout + "* Scan does not include images\n"
        
        return lout
    
    ################################################################
    def get_scaler(self,label=None):
        """
        return scaler
        """
        if label == None:
            label = self.primary_det[0]
        return self.scalers.get(label)
        
    ################################################################
    def get_positioner(self,label=None):
        """
        return positioner
        """
        if label == None:
            label = self.primary_axis[0]
        return self.positioners.get(label)

########################################################################
def append(data1,data2,sort=True):
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
def merge(data=[],average=True,align=False,fast=True):
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

