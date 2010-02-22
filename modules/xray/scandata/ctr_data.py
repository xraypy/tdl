"""
Functions for extracting ctr data from ScanData objects.

Authors / Modifications:
------------------------
T. Trainor (tptrainor@alaska.edu)
Frank Heberling (Frank.Heberling@ine.fzk.de)

Notes:
------



References:
-----------
 * E. Vlieg, J. Appl. Cryst. (1997). 30, 532-543
 * C. Schlepuetz et al, Acta Cryst. (2005). A61, 418-425

Todo:
----
 * Test!
 * Averaging/merging and merge statistics
 * Corrections for rocking scans
 * Add normalized F plots  - divide by |Fctr|,
    * need to pass delta_H, delta_K for non-rational surfaces...

"""
##############################################################################

import types, copy
import numpy as num
from matplotlib import pyplot

import plotter
from   mathutil import cosd, sind, tand
from   mathutil import arccosd, arcsind, arctand

import image_data
from   xtal.active_area import active_area
import gonio_psic 

##############################################################################
def ctr_data(scans,ctr=None,I=None,Inorm=None,Ierr=None,Ibgr=None,
             corr_params=None,scan_type=None):
    """
    create a ctr instance or add scan data to an existing instance
    """
    if ctr == None:
        # check for defaults
        if I==None: I ='I'
        if Inorm==None: Inorm='io'
        if Ierr==None: Inorm='Ierr',
        if Ibgr==None: Ibgr='Ibgr',
        if corr_params==None: corr_params={}
        if scan_type==None:scan_type='image'
        return CtrData(scans=scans,I=I,Inorm=Inorm,Ierr=Ierr,Ibgr=Ibgr,
                       corr_params=corr_params,scan_type=scan_type)
    else:
        ctr.append_scans(scans,I=I,Inorm=Inorm,Ierr=Ierr,Ibgr=Ibgr,
                         corr_params=corr_params,scan_type=scan_type)
    
##############################################################################
class CtrData:
    """
    CTR data

    scans = list of scan data instances

    I = string label corresponding to the intensity array, ie
        let y be the intensity, then y = scan[I]

    Inorm=string label corresponding to the normalization intensity array, ie
        norm = scan[Inorm], where the normalized intensity is taken as
           ynorm= y/norm

    Ierr = string label corresponding to the intensity error array, ie
        y_err = scan[Ierr].  We assume these are standard deviations 
        of the intensities (y).
        
      Note when the data are normalized we assume the statistics of norm 
        go as norm_err = sqrt(norm).  Therefore the error of normalized
        data is
           ynorm_err = ynorm * sqrt( (y_err/y)^2 + (norm_err/(norm))^2 )
                     = ynorm * sqrt( (y_err/y)^2 + 1/norm )
        If its assumed that y_err = sqrt(y) then the expression could be
        simplified futher, but we wont make that assumption since there
        may be other factors that determine how yerr was calculated.  

    Ibgr = string label corresponding to background intensity array

    corr_params = a dictionary containing the necessary information for
                  data corrections.
        corr_params['geom'] = Type of goniometer geometry ('psic' is default)
        corr_params['beam_slits'] = dictionary describing the beam slits,e.g.
                                    {'horz':.6,'vert':.8}
        corr_params['det_slits'] = dictionary describing the beam slits,e.g.
                                    {'horz':.6,'vert':.8}
        corr_params['sample'] = either dimater of round sample, or dictionary
                                describing the sample shape.
        corr_params['scale'] = scale factor to multiply by all the intensity
                               values. e.g.  if Io ~ 1million cps
                               then using 1e6 as the scale makes the normalized
                               intensity close to cps.  ie y = scale*y/norm

        See the Correction class for more info...

    scan_type = Type of scans (e.g. 'image', 'phi', etc..)
    
    """
    ##########################################################################
    def __init__(self,scans=[],I='I',Inorm='io',Ierr='Ierr',
                 Ibgr='Ibgr',corr_params={},scan_type='image'):
        """
        Constructor.
        
        """
        self.fig    = None
        self.cursor = None
        self.bad    = []
        self.scan   = []
        self.scan_index  = []
        #
        self.labels      = {'I':[],'Inorm':[],'Ierr':[],'Ibgr':[]}
        self.corr_params = []
        self.scan_type   = []
        #
        self.H     = num.array([],dtype=float)
        self.K     = num.array([],dtype=float)
        self.L     = num.array([],dtype=float)
        self.I     = num.array([],dtype=float)
        self.Inorm = num.array([],dtype=float)
        self.Ierr  = num.array([],dtype=float)
        self.Ibgr  = num.array([],dtype=float)
        self.ctot  = num.array([],dtype=float)
        self.F     = num.array([],dtype=float)
        self.Ferr  = num.array([],dtype=float)
        #
        self.append_scans(scans,I=I,Inorm=Inorm,Ierr=Ierr,Ibgr=Ibgr,
                          corr_params=corr_params,
                          scan_type=scan_type)

    ##########################################################################
    def __repr__(self,):
        """ """
        lout = "CTR DATA\n"
        lout = "%sNumber of scans = %i\n" % (lout,len(self.scan))
        lout = "%sNumber of structure factors = %i\n" % (lout,len(self.L))
        return lout

    ##########################################################################
    def __save__(self,):
        """ """
        del self.cursor
        self.cursor = None

    ##########################################################################
    def append_scans(self,scans,I=None,Inorm=None,Ierr=None,Ibgr=None,
                     corr_params=None,scan_type=None):
        """
        scans is a list of scan data objects.

        The rest of the arguments (defined above)
        should be the same for each scan in the list.

        For any argument with None passed we use previous vals...
        """
        if type(scans) != types.ListType:
            scans = [scans]
        
        # if None passed use the last values
        if I == None:           I = self.labels['I'][-1]
        if Inorm == None:       Inorm = self.labels['Inorm'][-1]
        if Ierr == None:        Ierr = self.labels['Ierr'][-1]
        if Ibgr == None:        Ibgr = self.labels['Ibgr'][-1]
        if corr_params == None: corr_params = self.corr_params[-1]
        if scan_type == None:   scan_type = self.scan_type[-1]

        # get all the data parsed out of each scan and append
        for scan in scans:
            data = self._scan_data(scan,I,Inorm,Ierr,Ibgr,corr_params,scan_type)
            if data == None: return

            #self.scan.append([])
            self.scan.append(scan)
            #
            self.scan_index.extend(data['scan_index'])
            self.labels['I'].extend(data['I_lbl'])
            self.labels['Inorm'].extend(data['Inorm_lbl'])
            self.labels['Ierr'].extend(data['Ierr_lbl'])
            self.labels['Ibgr'].extend(data['Ibgr_lbl'])
            self.corr_params.extend(data['corr_params'])
            self.scan_type.extend(data['scan_type'])
            #
            self.H     = num.append(self.H,data['H'])
            self.K     = num.append(self.K,data['K'])
            self.L     = num.append(self.L,data['L'])
            self.I     = num.append(self.I,data['I'])
            self.Inorm = num.append(self.Inorm,data['Inorm'])
            self.Ierr  = num.append(self.Ierr,data['Ierr'])
            self.Ibgr  = num.append(self.Ibgr,data['Ibgr'])
            self.ctot  = num.append(self.ctot,data['ctot'])
            self.F     = num.append(self.F,data['F'])
            self.Ferr  = num.append(self.Ferr,data['Ferr'])

    ##########################################################################
    def _scan_data(self,scan,I,Inorm,Ierr,Ibgr,corr_params,scan_type):
        """
        Parse scan into data...
        """
        data = {'scan_index':[],'I_lbl':[],'Inorm_lbl':[],
                'Ierr_lbl':[],'Ibgr_lbl':[],'corr_params':[],'scan_type':[],
                'H':[],'K':[],'L':[],'I':[],'Inorm':[],'Ierr':[],'Ibgr':[],
                'ctot':[],'F':[],'Ferr':[]}

        # compute a scan index
        scan_idx = len(self.scan)
        
        # image scan -> each scan point is a unique HKL
        if scan_type == 'image':
            if scan.image._is_integrated == False:
                scan.image.integrate()
            npts = int(scan.dims[0])
            for j in range(npts):
                data['scan_index'].append((scan_idx,j))
                data['I_lbl'].append(I)
                data['Inorm_lbl'].append(Inorm)
                data['Ierr_lbl'].append(Ierr)
                data['Ibgr_lbl'].append(Ibgr)
                data['corr_params'].append(corr_params)
                data['scan_type'].append(scan_type)
                #
                data['H'].append(scan['H'][j])
                data['K'].append(scan['K'][j])
                data['L'].append(scan['L'][j])
                # get F
                d = image_point_F(scan,j,I=I,Inorm=Inorm,
                                  Ierr=Ierr,Ibgr=Ibgr,
                                  corr_params=corr_params)
                data['I'].append(d['I'])
                data['Inorm'].append(d['Inorm'])
                data['Ierr'].append(d['Ierr'])
                data['Ibgr'].append(d['Ibgr'])
                data['ctot'].append(d['ctot'])
                data['F'].append(d['F'])
                data['Ferr'].append(d['Ferr'])
        return data

    ##########################################################################
    def integrate_point(self,idx,**kw):
        """
        (Re)-integrate an individual structure factor point.
        
        idx is the index number of the point

        If scan type is image use the following kw arguments
        (if the argument is not passed the existing value is used,
        ie just use these to update parameters)
        #
        bad        = True/False flags point as bad or not
        #
        roi        = image roi
        rotangle   = image rotation angle
        bgr_params = image background parameters
        plot       = True/False to show integration plot 
        fig        = Fig number for plot
        # 
        I          = Intensity label
        Inorm      = Intensity norm label
        Ierr       = Intensity error label
        Ibgr       = Intensity background label
        corr_params = CTR correction parameters
        
        """
        if idx not in range(len(self.L)): return

        bad = kw.get('bad')
        if bad != None:
            if bad == True:
                if idx not in self.bad:
                    self.bad.append(idx)
            elif bad == False:
                if idx in self.bad:
                    self.bad.remove(idx)
            else:
                print "Warning: bad should be True/False"

        if self.scan_type[idx]=="image":
            (scan_idx,point) = self.scan_index[idx]
            scan = self.scan[scan_idx]
            if scan.image._is_init() == False:
                scan.image._init_image()
            
            # parse integration parameters
            roi        = kw.get('roi')
            rotangle   = kw.get('rotangle')
            bgr_params = kw.get('bgr_params')
            #
            plot       = kw.get('plot',False)
            fig        = kw.get('fig')
            
            if idx in self.bad:
                bad = [point]
            else:
                bad = []
            # integrate the scan.  
            scan.image.integrate(idx=[point],roi=roi,rotangle=rotangle,
                                 bgr_params=bgr_params,bad_points=bad,
                                 plot=plot,fig=fig)
            
            # parse all the correction info and re-compute 
            I           = kw.get('I',self.labels['I'][idx])
            Inorm       = kw.get('Inorm',self.labels['Inorm'][idx])
            Ierr        = kw.get('Ierr', self.labels['Ierr'][idx])
            Ibgr        = kw.get('Ibgr', self.labels['Ibgr'][idx])
            corr_params = kw.get('corr_params',self.corr_params[idx])
            d = image_point_F(scan,point,I=I,Inorm=Inorm,
                              Ierr=Ierr,Ibgr=Ibgr,
                              corr_params=corr_params)
            
            # store results
            self.labels['I'][idx]     = I
            self.labels['Inorm'][idx] = Inorm
            self.labels['Ierr'][idx]  = Ierr
            self.corr_params[idx]     = corr_params
            self.H[idx]               = scan['H'][point]
            self.K[idx]               = scan['K'][point]
            self.L[idx]               = scan['L'][point]
            self.I[idx]               = d['I']
            self.Inorm[idx]           = d['Inorm']
            self.Ierr[idx]            = d['Ierr']
            self.Ibgr[idx]            = d['Ibgr']
            self.ctot[idx]            = d['ctot']
            self.F[idx]               = d['F']
            self.Ferr[idx]            = d['Ferr']
        return 

    ##########################################################################
    def plot(self,fig=None,num_col=2,cursor=True,verbose=True,spnt=None):
        """
        Plot the raw structure factor data
        """
        hksets  = sort_data(self)
        nset    = len(hksets)
        num_col = float(num_col)
        num_row = num.ceil(nset/num_col)
        pyplot.figure(fig)
        pyplot.clf()
        for j in range(nset):
            pyplot.subplot(num_row,num_col,j+1)
            d = hksets[j]
            title = 'H=%2.3f,K=%2.3f' % (d['H'][0],d['K'][0])
            pyplot.title(title, fontsize = 12)
            pyplot.plot(d['L'],d['F'],'b.-')
            pyplot.errorbar(d['L'],d['F'],d['Ferr'], fmt ='o')
            if spnt != None:
                if spnt in d['point_idx']:
                    idx = num.where(d['point_idx']==spnt)
                    pyplot.plot(d['L'][idx],d['F'][idx],'ro')
            pyplot.semilogy()
            #
            min_L = num.floor(num.min(d['L']))
            max_L = num.ceil(num.max(d['L']))
            idx   = num.where(d['F'] > 0.)
            min_F = num.min(d['F'][idx])
            min_F = 10.**(num.round(num.log10(min_F)) - 1)
            max_F = num.max(d['F'][idx])
            max_F = 10.**(num.round(num.log10(max_F)) + 1)
            pyplot.axis([min_L,max_L,min_F,max_F])
            #
            pyplot.xlabel('L')
            pyplot.ylabel('|F|')
        fig = pyplot.gcf()
        self.fig = fig.number
        if self.cursor != None:
            self.cursor._disconnect()
            del self.cursor
            self.cursor = None
        if cursor == True:
            self.cursor = plotter.cursor(fig=self.fig,verbose=verbose)

    ##########################################################################
    def plot_I(self,fig=None,num_col=2,cursor=True,verbose=True,spnt=None):
        """
        Plot the raw intensities
        """
        hksets  = sort_data(self)
        nset    = len(hksets)
        num_col = float(num_col)
        num_row = num.ceil(nset/num_col)
        pyplot.figure(fig)
        pyplot.clf()
        for j in range(nset):
            pyplot.subplot(num_row,num_col,j+1)
            d = hksets[j]
            title = 'H=%2.3f,K=%2.3f' % (d['H'][0],d['K'][0])
            pyplot.title(title, fontsize = 12)
            I  = d['I']
            In = d['Inorm']
            Ib = d['Ibgr']
            Ie = d['Ierr']
            y  = I/In
            yb = Ib/In
            ye = Ie/In
            #
            pyplot.plot(d['L'],y,'b.-',label='I/Inorm')
            pyplot.errorbar(d['L'],y,ye, fmt ='o')
            pyplot.plot(d['L'],yb,'m.-',label='Ibgr/Inorm')
            #
            min_L = num.floor(num.min(d['L']))
            max_L = num.ceil(num.max(d['L']))
            #
            idx   = num.where(y > 0.)
            min_I = num.min(y[idx])
            min_I = 10.**(num.round(num.log10(min_I)) - 1)
            idxb  = num.where(yb > 0.)
            if len(idxb[0]) == 0:
                min_Ibgr = min_I
            else:
                min_Ibgr = num.min(yb[idxb])
                min_Ibgr = 10.**(num.round(num.log10(min_Ibgr)) - 1)
            min_I = min(min_I,min_Ibgr)
            #
            max_I = num.max(y[idx])
            max_I = 10.**(num.round(num.log10(max_I)) + 1)
            if len(idxb[0]) == 0:
                max_Ibgr = max_I
            else:
                max_Ibgr = num.max(yb[idxb])
                max_Ibgr = 10.**(num.round(num.log10(max_Ibgr)) + 1)
            max_I = max(max_I,max_Ibgr)
            #
            idx  = num.where(I <= 0.)
            tmp  = I[idx] * 0.0 + 1.1*min_I 
            pyplot.plot(d['L'][idx],tmp,'bo')
            #
            if spnt != None:
                if spnt in d['point_idx']:
                    idx = num.where(d['point_idx']==spnt)
                    if I[idx] <= 0.:
                        pyplot.plot(d['L'][idx],[1.1*min_I],'ro')
                    else:
                        pyplot.plot(d['L'][idx],y[idx],'ro')
            #
            pyplot.semilogy()
            pyplot.axis([min_L,max_L,min_I,max_I])
            pyplot.xlabel('L')
            pyplot.ylabel('Intensity')
            pyplot.legend()
        fig = pyplot.gcf()
        self.fig = fig.number
        if self.cursor != None:
            self.cursor._disconnect()
            del self.cursor
            self.cursor = None
        if cursor == True:
            self.cursor = plotter.cursor(fig=self.fig,verbose=verbose)

    ##########################################################################
    def plot_point(self,idx=None,fig=None,show_int=False,cmap=None):
        """
        Plot the raw data for a selected point

        idx = point index.
        if idx = None, then uses last cursor click
        fig = fig number
        """
        if idx == None:
            idx = self.get_idx()
        if self.scan_type[idx] == 'image':
            if show_int:
                self.integrate_point(idx,plot=True,fig=fig)
            else:
                (scan_idx,point) = self.scan_index[idx]
                self.scan[scan_idx].image.plot(idx=point,fig=fig,cmap=cmap)
        else:
            # plot scan data
            pass

    ##########################################################################
    def get_idx(self):
        """
        Get point index from plot selection
        
        ie  L = ctr.L[point_idx]
        etc.
        """
        if self.cursor == None:
            return None
        if self.cursor.clicked == False:
            #self.cursor.get_click(msg="Select a data point")
            return None
        L = self.cursor.x
        subplot = self.cursor.subplot
        if subplot < 0:
            return None
        hksets  = sort_data(self)
        return self._get_idx(subplot,L,hksets)
    
    def _get_idx(self,subplot,L,hksets):
        d   = hksets[subplot]
        tmp = num.fabs(d['L']-L)
        idx = num.where(tmp==min(tmp))
        if len(idx) > 0:
            idx = idx[0]
        point_idx = d['point_idx'][idx]
        return point_idx

    ##########################################################################
    def get_scan(self,idx=None):
        """
        Get scan from point index

        Returns (scan,point)
        """
        if idx == None:
            idx = self.get_idx()
            if idx == None: return None
        scan_idx = self.scan_index[idx]
        if len(scan_idx) > 1:
            point = scan_idx[1]
            idx = scan_idx[0]
        else:
            idx = scan_idx
            point = 0
        scan = self.scan[idx]
        return (scan,point)
    
    ##########################################################################
    def get_points(self,fig=None):
        """
        get index value of points from plot
        """
        if self.cursor == None:
            return None
        if self.cursor.clicked == False:
            #self.cursor.get_click(msg="Zoom on plot")
            return None
        self.cursor._zoom()
        z = self.cursor.zoom
        z = (z[0][0],z[1][0])
        Lmin = min(z)
        Lmax = max(z)
        subplot = self.cursor.subplot
        if subplot < 0:
            return None
        hksets  = sort_data(self)
        return self._get_idx_range(subplot,Lmin,Lmax,hksets)
        
    def _get_idx_range(self,subplot,Lmin,Lmax,hksets):
        d    = hksets[subplot]
        tmp  = d['L']
        idx1 = num.where(tmp<Lmin)
        idx2 = num.where(tmp>Lmax)
        tmp  = num.fabs(tmp)
        tmp[idx1] = 0.0
        tmp[idx2] = 0.0
        idx = num.where(tmp > 0 )
        point_idx = d['point_idx'][idx]
        return (point_idx)
    
    ##########################################################################
    def write_HKL(self,fname = 'ctr.lst'):
        """
        dump data file
        """
        f = open(fname, 'w')
        header = "#idx %5s %5s %5s %7s %7s\n" % ('H','K','L','F','Ferr')
        f.write(header)
        for i in range(len(self.L)):
            if self.I[i] > 0:
                line = "%4i %3.2f %3.2f %6.3f %6.6g %6.6g\n" % (i,round(self.H[i]),
                                                                round(self.K[i]),
                                                                self.L[i],self.F[i],
                                                                self.Ferr[i])
                f.write(line)
        f.close()

##########################################################################
def sort_data(ctr,hkdecimal=3):
    """
    Return a dict of sorted data

    Assume H,K define a set with a range of L values
    All arrays should be of len npts. 

    """
    # round H and K to sepcified precision
    H = num.around(ctr.H,decimals=hkdecimal)
    K = num.around(ctr.K,decimals=hkdecimal)
    L = ctr.L
    F = ctr.F
    Ferr = ctr.Ferr
    #
    I     = ctr.I
    Inorm = ctr.Inorm
    Ierr  = ctr.Ierr
    Ibgr  = ctr.Ibgr
    #
    scan_idx = ctr.scan_index
    npts = len(F)

    #find all unique sets
    hkvals = []
    for j in range(npts):
        s = (H[j],K[j]) 
        if s not in hkvals:
            hkvals.append(s)

    # sort the hkvals
    # and stick data in correct set
    hkvals.sort()
    nsets = len(hkvals)
    #d = {'H':[],'K':[],'L':[],'F':[],'Ferr':[],'idx':[]}
    #hkset  = [copy.copy(d) for j in range(nsets)]
    hkset = []
    for j in range(nsets):
        hkset.append({'H':[],'K':[],'L':[],'F':[],'Ferr':[],
                      'I':[],'Inorm':[],'Ierr':[],'Ibgr':[],
                      'point_idx':[],'scan_idx':[]})

    for j in range(npts):
        s      = (H[j],K[j])
        setidx = hkvals.index(s)
        hkset[setidx]['H'].append(H[j])
        hkset[setidx]['K'].append(K[j])
        hkset[setidx]['L'].append(L[j])
        hkset[setidx]['F'].append(F[j])
        hkset[setidx]['Ferr'].append(Ferr[j])
        hkset[setidx]['I'].append(I[j])
        hkset[setidx]['Inorm'].append(Inorm[j])
        hkset[setidx]['Ierr'].append(Ierr[j])
        hkset[setidx]['Ibgr'].append(Ibgr[j])
        hkset[setidx]['point_idx'].append(j)
        hkset[setidx]['scan_idx'].append(scan_idx[j])

    # make arrays num arrays
    for j in range(nsets):
        hkset[j]['H'] = num.array(hkset[j]['H'])
        hkset[j]['K'] = num.array(hkset[j]['K'])
        hkset[j]['L'] = num.array(hkset[j]['L'])
        hkset[j]['F'] = num.array(hkset[j]['F'])
        hkset[j]['Ferr'] = num.array(hkset[j]['Ferr'])
        hkset[j]['I']     = num.array(hkset[j]['I'])
        hkset[j]['Inorm'] = num.array(hkset[j]['Inorm'])
        hkset[j]['Ierr']  = num.array(hkset[j]['Ierr'])
        hkset[j]['Ibgr']  = num.array(hkset[j]['Ibgr'])
        hkset[j]['point_idx']  = num.array(hkset[j]['point_idx'])
        hkset[j]['scan_idx']  = num.array(hkset[j]['scan_idx'])

    # now sort each set by L
    for j in range(nsets):
        lidx = num.argsort(hkset[j]['L'])
        hkset[j]['H'] = hkset[j]['H'][lidx]
        hkset[j]['K'] = hkset[j]['K'][lidx]
        hkset[j]['L'] = hkset[j]['L'][lidx]
        hkset[j]['F'] = hkset[j]['F'][lidx]
        hkset[j]['Ferr'] = hkset[j]['Ferr'][lidx]
        hkset[j]['I'] = hkset[j]['I'][lidx]
        hkset[j]['Inorm'] = hkset[j]['Inorm'][lidx]
        hkset[j]['Ierr'] = hkset[j]['Ierr'][lidx]
        hkset[j]['Ibgr'] = hkset[j]['Ibgr'][lidx]
        hkset[j]['point_idx'] = hkset[j]['point_idx'][lidx]
        hkset[j]['scan_idx'] = hkset[j]['scan_idx'][lidx]

    return hkset

##############################################################################
def image_point_F(scan,point,I='I',Inorm='io',Ierr='Ierr',Ibgr='Ibgr',
                  corr_params={}):
    """
    compute F for a single scan point in an image scan
    """
    d = {'I':0.0,'Inorm':0.0,'Ierr':0.0,'Ibgr':0.0,'ctot':1.0,'F':0.0,'Ferr':0.0}
    d['I']     = scan[I][point]
    d['Inorm'] = scan[Inorm][point]
    d['Ierr']  = scan[Ierr][point]
    d['Ibgr']  = scan[Ibgr][point]
    if corr_params == None:
        d['ctot'] = 1.0
        scale = 1.0
    else:
        # compute correction factors
        scale  = corr_params.get('scale',1.0)
        scale  = float(scale)
        corr = _get_corr(scan,point,corr_params)
        if corr == None:
            d['ctot'] = 1.0
        else:
            d['ctot'] = corr.ctot_stationary()
    # compute F
    if d['I'] <= 0.0 or d['Inorm'] <= 0.:
        d['F']    = 0.0
        d['Ferr'] = 0.0
    else:
        yn     = scale*d['I']/d['Inorm']
        yn_err = yn * num.sqrt( (d['Ierr']/d['I'])**2. + 1./d['Inorm'] )
        d['F']    = num.sqrt(d['ctot']*yn)
        d['Ferr'] = num.sqrt(d['ctot']*yn_err)
    
    return d

##############################################################################
def _get_corr(scan,point,corr_params):
    """
    get CtrCorrection instance
    """
    geom   = corr_params.get('geom')
    if geom == None: geom='psic'
    beam   = corr_params.get('beam_slits',{})
    det    = corr_params.get('det_slits')
    sample = corr_params.get('sample')
    # get gonio instance for corrections
    if geom == 'psic':
        gonio = gonio_psic.psic_from_spec(scan['G'])
        _update_psic_angles(gonio,scan,point)
        corr  = CtrCorrectionPsic(gonio=gonio,beam_slits=beam,
                                  det_slits=det,sample=sample)
    else:
        print "Geometry %s not implemented" % geom
        corr = None
    return corr

##############################################################################
def get_params(ctr,point):
    """
    Return relevant parameters from a specified point
    of a ctr object.  ie use to copy parameters...
    
    """
    intpar = {}
    corrpar = {}
    (scan,spnt) = ctr.get_scan(point)
    #
    intpar['I']     = ctr.labels['I'][point]
    intpar['Inorm'] = ctr.labels['Inorm'][point]
    intpar['Ierr']  = ctr.labels['Ierr'][point]
    intpar['Ibgr']  = ctr.labels['Ibgr'][point]
    if ctr.scan_type[point] == 'image':
        intpar['image roi']      = scan.image.rois[spnt]
        intpar['image rotangle'] = scan.image.rotangle[spnt]
        intpar['bgr flag']       = scan.image.bgrpar[spnt]['bgrflag']
        intpar['bgr col nbgr']   = scan.image.bgrpar[spnt]['cnbgr']
        intpar['bgr col width']  = scan.image.bgrpar[spnt]['cwidth']
        intpar['bgr col power']  = scan.image.bgrpar[spnt]['cpow']
        intpar['bgr col tan']    = scan.image.bgrpar[spnt]['ctan']
        intpar['bgr row nbgr']   = scan.image.bgrpar[spnt]['rnbgr']
        intpar['bgr row width']  = scan.image.bgrpar[spnt]['rwidth']
        intpar['bgr row power']  = scan.image.bgrpar[spnt]['rpow']
        intpar['bgr row tan']    = scan.image.bgrpar[spnt]['rtan']
    else:
        pass
    #
    corrpar['beam_slits'] = ctr.corr_params[point].get('beam_slits')
    corrpar['det_slits'] = ctr.corr_params[point].get('det_slits')
    corrpar['geom'] = ctr.corr_params[point].get('geom')
    corrpar['scale'] = ctr.corr_params[point].get('scale')
    sample  = ctr.corr_params[point].get('sample')
    if sample == None:
        corrpar['sample dia']     = None
        corrpar['sample polygon'] = None
        corrpar['sample angles']  = None
    elif type(sample) == types.DictType:
        corrpar['sample dia']     = None
        corrpar['sample polygon'] = ctr.corr_params[point]['sample']['polygon']
        corrpar['sample angles']  = ctr.corr_params[point]['sample']['angles']
    else:
        corrpar['sample dia']     = ctr.corr_params[point]['sample']
        corrpar['sample polygon'] = None
        corrpar['sample angles']  = None

    return (intpar,corrpar)

##############################################################################
def set_params(ctr,point,intpar={},corrpar={}):
    """
    Set ctr parameters.

    The intpar and corrpar arguments should be from get_param fcn above
    ie that sets the correct format...
    """
    (scan,spnt) = ctr.get_scan(point)
    #
    if len(intpar) > 0:
        ctr.labels['I'][point]     = intpar['I']
        ctr.labels['Inorm'][point] = intpar['Inorm']
        ctr.labels['Ierr'][point]  = intpar['Ierr']
        ctr.labels['Ibgr'][point]  = intpar['Ibgr']
        if ctr.scan_type[point] == 'image':
            if type(intpar['image roi']) == types.StringType:
                scan.image.rois[spnt] = eval(intpar['image roi'])
            else:
                scan.image.rois[spnt] = intpar['image roi']
            if type(intpar['image rotangle']) == types.StringType:
                scan.image.rotangle[spnt] = eval(intpar['image rotangle'])
            else:
                scan.image.rotangle[spnt] = intpar['image rotangle']
            if type(intpar['bgr flag']) == types.StringType:
                scan.image.bgrpar[spnt]['bgrflag'] = eval(intpar['bgr flag'])
            else:
                scan.image.bgrpar[spnt]['bgrflag'] = intpar['bgr flag']
            if type(intpar['bgr col nbgr']) == types.StringType:
                scan.image.bgrpar[spnt]['cnbgr'] = eval(intpar['bgr col nbgr'])
            else:
                scan.image.bgrpar[spnt]['cnbgr'] = intpar['bgr col nbgr']
            if type(intpar['bgr col width']) == types.StringType:
                scan.image.bgrpar[spnt]['cwidth'] = eval(intpar['bgr col width'])
            else:
                scan.image.bgrpar[spnt]['cwidth'] = intpar['bgr col width']
            if type(intpar['bgr col power']) == types.StringType:
                scan.image.bgrpar[spnt]['cpow'] = eval(intpar['bgr col power'])
            else:
                scan.image.bgrpar[spnt]['cpow'] = intpar['bgr col power']
            if type(intpar['bgr col tan']) == types.StringType:
                scan.image.bgrpar[spnt]['ctan'] = eval(intpar['bgr col tan'])
            else:
                scan.image.bgrpar[spnt]['ctan'] = intpar['bgr col tan']
            if type(intpar['bgr row nbgr']) == types.StringType:
                scan.image.bgrpar[spnt]['rnbgr'] = eval(intpar['bgr row nbgr'])
            else:
                scan.image.bgrpar[spnt]['rnbgr'] = intpar['bgr row nbgr']
            if type(intpar['bgr row width']) == types.StringType:
                scan.image.bgrpar[spnt]['rwidth'] = eval(intpar['bgr row width'])
            else:
                scan.image.bgrpar[spnt]['rwidth'] = intpar['bgr row width']
            if type(intpar['bgr row power']) == types.StringType:
                scan.image.bgrpar[spnt]['rpow'] = eval(intpar['bgr row power'])
            else:
                scan.image.bgrpar[spnt]['rpow'] = intpar['bgr row power']
            if type(intpar['bgr row tan']) == types.StringType:
                scan.image.bgrpar[spnt]['rtan'] = eval(intpar['bgr row tan'])
            else:
                scan.image.bgrpar[spnt]['rtan'] = intpar['bgr row tan']
    #
    if len(corrpar) > 0:
        ctr.corr_params[point] = {}
        for (key,val) in corrpar.items():
            if type(corrpar['beam_slits']) == types.StringType:
                ctr.corr_params[point]['beam_slits'] = eval(corrpar['beam_slits'])
            else:
                ctr.corr_params[point]['beam_slits'] = corrpar['beam_slits']
            if type(corrpar['det_slits']) == types.StringType:
                ctr.corr_params[point]['det_slits'] = eval(corrpar['det_slits'])
            else:
                ctr.corr_params[point]['det_slits'] = corrpar['det_slits']
            ctr.corr_params[point]['geom'] = corrpar['geom']
            if type(corrpar['scale']) == types.StringType:
                ctr.corr_params[point]['scale'] = eval(corrpar['scale'])
            else:
                ctr.corr_params[point]['scale'] = corrpar['scale']
            #
            if type(corrpar['sample dia']) == types.StringType:
                sdia = eval(corrpar['sample dia'])
            else:
                sdia = corrpar['sample dia']
            if type(corrpar['sample polygon']) == types.StringType:
                spoly = eval(corrpar['sample polygon'])
            else:
                spoly = corrpar['sample polygon']
            if type(corrpar['sample angles']) == types.StringType:
                sangles = eval(corrpar['sample angles'])
            else:
                sangles = corrpar['sample angles']
            if sdia != None:
                ctr.corr_params[point]['sample'] = sdia
            elif spoly != None:
                ctr.corr_params[point]['sample'] = {}
                ctr.corr_params[point]['sample']['polygon'] = spoly
                ctr.corr_params[point]['sample']['angles'] = sangles
            else:
                ctr.corr_params[point]['sample'] = None
    
##############################################################################
def _update_psic_angles(gonio,scan,point):
    """
    given a psic gonio instance, a scandata object
    and a scan point, update the gonio angles...
    """
    npts = int(scan.dims[0])
    if len(scan['phi']) == npts:
        phi=scan['phi'][point]
    else:
        phi=scan['phi']
    if len(scan['chi']) == npts:
        chi=scan['chi'][point]
    else:
        chi=scan['chi']
    if len(scan['eta']) == npts:
        eta=scan['eta'][point]
    else:
        eta=scan['eta']
    if len(scan['mu']) == npts:
        mu=scan['mu'][point]
    else:
        mu=scan['mu']
    if len(scan['nu']) == npts:
        nu=scan['nu'][point]
    else:
        nu=scan['nu']
    if len(scan['del']) == npts:
        delta=scan['del'][point]
    else:
        delta=scan['del']
    #
    gonio.set_angles(phi=phi,chi=chi,eta=eta,
                     mu=mu,nu=nu,delta=delta)

##############################################################################
class CtrCorrectionPsic:
    """
    Data point operations / corrections for Psic geometry

    Note: All correction factors are defined such that the
    measured data is corrected by multiplying times
    the correction: 
      Ic  = Im*ct
    where
      Im = Idet/Io = uncorrected (measured) intensity

    In other words we use the following formalism:
      Im = (|F|**2)* prod_i(Xi)
    where Xi are various (geometric) factors that 
    influence the measured intensity.  To get the
    structure factor:
      |F| = sqrt(Im/prod_i(Xi)) = sqrt(Im* ct)
    and
      ct = prod_i(1/Xi) = prod_i(ci)
      ci = 1/Xi
      
    If there is an error or problem in the routine for a sepcific
    correction factor, (e.g. divide by zero), the routine should
    return a zero.  This way the corrected data is zero'd....

    The correction factors depend on the goniometer geometry
     gonio = gonio_psic.Psic instance

    The slits settings are needed.  Note if using a large area detector
    you may pass det_slits = None and just spill off will be computed
       beam_slits = {'horz':.6,'vert':.8}
       det_slits = {'horz':20.0,'vert':10.5}
    these are defined wrt psic phi-frame:
       horz = beam/detector horz width (total slit width in lab-z,
              or the horizontal scattering plane)
       vert = detector vert hieght (total slit width in lab-x,
              or the vertical scattering plane)

    A sample description is needed.
    If sample = a number then is is taken as the diameter
    of a round sample mounted on center.

    Otherwise is may describe a general sample shape:    
      sample = {}
        sample['polygon'] = [[1.,1.], [.5,1.5], [-1.,1.],
                             [-1.,-1.],[0.,.5],[1.,-1.]]
        sample['angles']  = {'phi':108.0007,'chi':0.4831}

        polygon = [[x,y,z],[x,y,z],[x,y,z],....]
                  is a list of vectors that describe the shape of
                  the sample.  They should be given in general lab
                  frame coordinates.

         angles = {'phi':0.,'chi':0.,'eta':0.,'mu':0.}
             are the instrument angles at which the sample
             vectors were determined.

      The lab frame coordinate systems is defined such that:
        x is vertical (perpendicular, pointing to the ceiling of the hutch)
        y is directed along the incident beam path
        z make the system right handed and lies in the horizontal scattering plane
          (i.e. z is parallel to the phi axis)

        The center (0,0,0) of the lab frame is the rotation center of the instrument.

        If the sample vectors are given at the flat phi and chi values and with
        the correct sample hieght (sample Z set so the sample surface is on the
        rotation center), then the z values of the sample vectors will be zero.
        If 2D vectors are passed we therefore assume these are [x,y,0].  If this
        is the case then make sure:
            angles = {'phi':flatphi,'chi':flatchi,'eta':0.,'mu':0.}

        The easiest way to determine the sample coordinate vectors is to take a picture
        of the sample with a camera mounted such that is looks directly down the omega
        axis and the gonio angles set at the sample flat phi and chi values and
        eta = mu = 0. Then find the sample rotation center and measure the position
        of each corner (in mm) with up being the +x direction, and downstream
        being the +y direction.  

    Note need to add attenuation parameters.  
    
    """
    def __init__(self,gonio=None,beam_slits={},det_slits=None,sample=None):
        self.gonio      = gonio
        if self.gonio.calc_psuedo == False:
            self.gonio.calc_psuedo = True
            self.gonio._update_psuedo()
        self.beam_slits = beam_slits
        self.det_slits  = det_slits
        self.sample     = sample
        # fraction horz polarization
        self.fh         = 1.0

    ##########################################################################
    def ctot_stationary(self,plot=False,fig=None):
        """
        correction factors for stationary measurements (e.g. images)
        """
        cp = self.polarization()
        cl = self.lorentz_stationary()
        ca = self.active_area(plot=plot,fig=fig)
        ct = (cp)*(cl)*(ca)
        if plot == True:
            print "Correction factors (mult by I)" 
            print "   Polarization=%f" % cp
            print "   Lorentz=%f" % cl
            print "   Area=%f" % ca
            print "   Total=%f" % ct
        return ct

    ##########################################################################
    def lorentz_stationary(self):
        """
        Compute the Lorentz factor for a stationary (image)
        measurement.  See Vlieg 1997

        Measured data is corrected for Lorentz factor as: 
          Ic  = Im * cl
        """
        beta  = self.gonio.pangles['beta']
        cl = sind(beta)
        return cl

    ##########################################################################
    def lorentz_scan(self):
        """
        Compute the Lorentz factor for a generic scan
        Note this is approximate. In general for bulk xrd
        with single crystals lorentz factor is defined as:
          L = 1/sin(2theta)
        We need L for specific scan types, e.g. phi, omega, etc..
        See Vlieg 1997

        Measured data is corrected for Lorentz factor as: 
          Ic  = Im * cl = Im/L
        """
        tth  = self.gonio.pangles['tth']
        cl = sind(tth)
        return cl

    ##########################################################################
    def rod_intercept(self,):
        """
        Compute the dl of the rod intercepted by the detector.
        (this correction only applies for rocking scans)
        This can be (roughly) approximated from:
          dl = dl_o * cos(beta)
          where dl_o is the detector acceptance at beta = 0

        Note this is an approximation for all but simple specular scans,
        Should use the detector acceptance polygon and determine
        the range of dl for the specific scan axis used.
        
        Measured data is corrected for the intercept as: 
          Ic  = Im * cr = Im/dl
        """
        beta  = self.gonio.pangles['beta']
        cr   = cosd(beta)
        return cr
    
    ##########################################################################
    def polarization(self,):
        """
        Compute polarization correction factor.
        
        For a horizontally polarized beam (polarization vector
        parrallel to the lab-frame z direction) the polarization
        factor is normally defined as:
           p = 1-(cos(del)*sin(nu))^2
        For a beam with mixed horizontal and vertical polarization:
           p = fh( 1-(cos(del)*sin(nu))^2 ) + (1-fh)(1-sin(del)^2)
        where fh is the fraction of horizontal polarization.

        Measured data is corrected for polarization as: 
          Ic  = Im * cp = Im/p
        """
        fh    = self.fh
        delta = self.gonio.angles['delta']
        nu    = self.gonio.angles['nu']
        p = 1. - ( cosd(delta) * sind(nu) )**2.
        if fh != 1.0:
            p = fh * c_p + (1.-fh)*(1.0 - (sind(delta))**2.)
        if p == 0.:
            cp = 0.
        else:
            cp = 1./p

        return cp

    ##########################################################################
    def active_area(self,plot=False,fig=None):
        """
        Compute active area correction (c_a = A_beam/A_int)
        Use to correct scattering data for area effects,
        including spilloff, i.e.
           Ic = Im * ca = Im/A_ratio 
           A_ratio = A_int/A_beam 
        where
           A_int = intersection area (area of beam on sample
                   viewed by detector)
           A_beam = total beam area
        """
        if self.beam_slits == {} or self.beam_slits == None:
            print "Warning beam slits not specified"
            return 1.0
        alpha = self.gonio.pangles['alpha']
        beta  = self.gonio.pangles['beta']
        if plot == True:
            print 'Alpha = ', alpha, ', Beta = ', beta
        if alpha < 0.0:
            print 'alpha is less than 0.0'
            return 0.0
        elif beta < 0.0:
            print 'beta is less than 0.0'
            return 0.0

        # get beam vectors
        bh = self.beam_slits['horz']
        bv = self.beam_slits['vert']
        beam = gonio_psic.beam_vectors(h=bh,v=bv)

        # get det vectors
        if self.det_slits == None:
            det = None
        else:
            dh = self.det_slits['horz']
            dv = self.det_slits['vert']
            det  = gonio_psic.det_vectors(h=dh,v=dv,
                                          nu=self.gonio.angles['nu'],
                                          delta=self.gonio.angles['delta'])
        # get sample poly
        if type(self.sample) == types.DictType:
            sample_vecs = self.sample['polygon']
            sample_angles = self.sample['angles']
            sample = gonio_psic.sample_vectors(sample_vecs,
                                               angles=sample_angles,
                                               gonio=self.gonio)
        else:
            sample = self.sample

        # compute active_area
        (A_beam,A_int) = active_area(self.gonio.nm,ki=self.gonio.ki,
                                     kr=self.gonio.kr,beam=beam,det=det,
                                     sample=sample,plot=plot,fig=fig)
        if A_int == 0.:
            ca = 0.
        else:
            ca = A_beam/A_int

        return ca

##############################################################################
##############################################################################
def test1():
    #psic = psic_from_spec(G,angles={})
    psic = gonio_psic.test2(show=False)
    psic.set_angles(phi=12.,chi=30.,eta=20.,
                    mu=25.,nu=75.,delta=20.)
    beam_slits = {'horz':.6,'vert':.8}
    det_slits = {'horz':20.0,'vert':10.5}
    sample = {}
    sample['polygon'] = [[1.,1.], [.5,1.5], [-1.,1.], [-1.,-1.],[0.,.5],[1.,-1.]]
    sample['angles']  = {'phi':108.0007,'chi':0.4831}
    cor = CtrCorrectionPsic(gonio=psic,beam_slits=beam_slits,
                            det_slits=det_slits,sample=sample)
    ct = cor.ctot_stationary(plot=True)

##########################################################################
if __name__ == "__main__":
    """
    test 
    """
    test1()



