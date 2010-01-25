##############################################################################
"""
T. Trainor (tptrainor@alaska.edu)
Frank Heberling (Frank.Heberling@ine.fzk.de)

Functions for extracting ctr data from ScanData objects

Modifications:
--------------
"""
##############################################################################
"""
Notes
------

References:
Vlieg 1997
Schlepuetz 2005


Todo
----
- Test!
- sorting and HK panel plots
- averaging/merging and merge statistics
- corrections for rocking scans
- editing / re-integrating select points:
  - select from plot and give type (image/scan)
"""
##############################################################################

import types
import numpy as num
from matplotlib import pyplot

from mathutil import cosd, sind, tand
from mathutil import arccosd, arcsind, arctand

from xtal.active_area import active_area
import gonio_psic 

##############################################################################
def ctr_data(scans,ctr=None,I=None,Inorm=None,Ierr=None,
             corr_params=None,scan_type=None):
    """
    creat a ctr instance or add scan data to an existing instance
    """
    if ctr == None:
        # check for defaults
        if I==None: I ='I_c'
        if Inorm==None: Inorm='Io'
        if Ierr==None: Inorm='Ierr_c',
        if corr_params==None: corr_params=[]
        if scan_type==None:scan_type='image'
        return CtrData(scans=scans,I=I,Inorm=Inorm,Ierr=Ierr,
                       corr_params=corr_params,scan_type=scan_type)
    else:
        ctr.append_data(scans,I=I,Inorm=Inorm,Ierr=Ierr,
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

    corr_params = a dictionary containing the necessary information for
                  data corrections.
        corr_params['geom'] = Type of goniometer geometry ('psic' is default)
        corr_params['beam_slits'] = dictionary describing the beam slits,e.g.
                                    {'horz':.6,'vert':.8}
        corr_params['det_slits'] = dictionary describing the beam slits,e.g.
                                    {'horz':.6,'vert':.8}
        corr_params['sample'] = either dimater of round sample, or dictionary
                                describing the sample shape.
        See the Correction class for more info...

    scan_type = Type of scans (e.g. 'image', 'phi', etc..)
    
    """
    def __init__(self,scans=[],I='I_c',Inorm='Io',Ierr='Ierr_c',
                 corr_params=[],scan_type='image'):
        #
        self.scan_data   = []
        self.scan_index  = []
        #
        self.labels      = {'I':[],'Inorm':[],'Ierr':[]}
        self.corr_params = []
        self.scan_type   = []
        #
        self.H     = num.array([],dtype=float)
        self.K     = num.array([],dtype=float)
        self.L     = num.array([],dtype=float)
        self.I     = num.array([],dtype=float)
        self.Inorm = num.array([],dtype=float)
        self.Ierr  = num.array([],dtype=float)
        self.ctot  = num.array([],dtype=float)
        self.F     = num.array([],dtype=float)
        self.Ferr  = num.array([],dtype=float)
        #
        self.append_data(scans,I=I,Inorm=Inorm,Ierr=Ierr,
                         corr_params=corr_params,
                         scan_type=scan_type)

    ##########################################################################
    def append_data(self,scans,I=None,Inorm=None,Ierr=None,
                    corr_params=None,scan_type=None):
        """
        scans is a list of scan data objects.  The rest of the arguments
        should be the same for each scan in the list or we use previous vals...
        """
        if type(scans) != types.ListType:
            scans = [scans]
        
        # if None passed use the last values
        if I == None: I = self.labels['I'][-1]
        if Inorm == None: Inorm = self.labels['Inorm'][-1]
        if Ierr == None: Ierr = self.labels['Ierr'][-1]
        if corr_params == None: corr_params = self.corr_params[-1]
        if scan_type == None: scan_type = self.scan_type[-1]

        # get all the data parsed out of each scan and append
        for scan in scans:
            data = _get_data(scan,I,Inorm,Ierr,corr_params,scan_type)
            if data == None: return

            # note we need to put images into hdf file to store..
            #self.scan_data.append(scan)
            self.scan_data.append([])
            #
            self.scan_index.extend(data['scan_index'])
            self.labels['I'].extend(data['I_lbl'])
            self.labels['Inorm'].extend(data['Inorm_lbl'])
            self.labels['Ierr'].extend(data['Ierr_lbl'])
            self.corr_params.extend(data['corr_params'])
            self.scan_type.extend(data['scan_type'])
            #
            self.H     = num.append(self.H,data['H'])
            self.K     = num.append(self.K,data['K'])
            self.L     = num.append(self.L,data['L'])
            self.I     = num.append(self.I,data['I'])
            self.Inorm = num.append(self.Inorm,data['Inorm'])
            self.Ierr  = num.append(self.Ierr,data['Ierr'])
            self.ctot  = num.append(self.ctot,data['ctot'])
            self.F     = num.append(self.F,data['F'])
            self.Ferr  = num.append(self.Ferr,data['Ferr'])

    ##########################################################################
    def _get_data(self,scan,I_lbl,Ierr_lbl,Inorm_lbl,corr_params,scan_type):
        """
        parse scan into data...
        """
        data = {'scan_index':[],'I_lbl':[],'Inorm_lbl':[],
                'Ierr_lbl':[],'corr_params':[],'scan_type':[],
                'H':[],'K':[],'L':[],'I':[],'Inorm':[],'Ierr':[],
                'ctot':[],'F':[],'Ferr':[]}

        # compute a scan index
        scan_idx = len(self.scan_data)
        
        # image scan -> each scan point is a unique HKL
        if scan_type == 'image':
            npts = int(scan.dims[0])
            for j in range(npts):
                data['scan_index'].append((scan_idx,j))
                data['I_lbl'].append(I_lbl)
                data['Inorm_lbl'].append(Inorm_lbl)
                data['Ierr_lbl'].append(Ierr_lbl)
                data['corr_params'].append(corr_params)
                data['scan_type'].append(scan_type)
                #
                data['H'].append(scan['H'][j])
                data['K'].append(scan['K'][j])
                data['L'].append(scan['L'][j])
                # get F
                (I,Inorm,Ierr,ctot,F,Ferr) = _image_point_F(scan,point,I_lbl,
                                                            Inorm_lbl,Ierr_lbl,
                                                            corr_params)
                data['I'].append(I)
                data['Inorm'].append(Inorm)
                data['Ierr'].append(Ierr)
                data['ctot'].append(ctot)
                data['F'].append(F)
                data['Ferr'].append(Ferr)
    
    ##########################################################################
    def _image_point_F(self,scan,point,I,Inorm,Ierr,corr_params):
        """
        compute F for a single scan point in an image scan
        """
        y     = scan[I][point]
        norm  = scan[Inorm][point]
        y_err = scan[Ierr][point]
        if corr_params == None:
            ctot = 1.0
        else:
            # compute correction factors
            geom   = corr_params.get('geom')
            if geom == None: geom='psic'
            beam   = corr_params.get('beam_slits')
            det    = corr_params.get('det_slits')
            sample = corr_params.get('sample')
            # get gonio instance for corrections
            if geom == 'psic':
                gonio = gonio_psic.psic_from_spec(scan['G'])
                update_psic_angles(gonio,scan,point)
                corr  = CtrCorrectionPsic(gonio=gonio,beam_slits=beam,
                                          det_slits=det,sample=sample)
            # else another geom or an error...
            ctot  = corr.ctot_stationary()
        # compute F
        yn     = y/norm
        yn_err = yn * num.sqrt( (y_err/y)**2. + 1./norm )

        F    = num.sqrt(ctot*yn)
        Ferr = num.sqrt(ctot*yn_err)
        
        return (y,norm,yn_err,ctot,F,Ferr)

    ##########################################################################
    def plot(self):
        """
        plot the rod.
        Need multi panel plot 
        """
        pyplot.figure()
        pyplot.plot(self.L, self.F,'b')
        pyplot.errorbar(self.L,self.F,self.Ferr, fmt ='o')
        pyplot.semilogy()
        pyplot.show()

    ##########################################################################
    def write_HKL(self,fname = 'ctr.hkl'):
        #write data file for row data
        f = open(fname, 'w')
        f.write('#idx        H            K            L            F           F_err \n')
        for i in range(len(self.L)):
            if self.I[i] > 0:
                line = "%i    %i    %i    %6.3f    %6.6g    %6.6g\n" % (i,
                                                                  round(self.H[i]),
                                                                  round(self.K[i]),
                                                                  self.L[i],
                                                                  self.F[i],
                                                                  self.Ferr[i])
                f.write(line)
        f.close()

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
    def __init__(self,gonio=None,beam_slits={},det_slits={},sample={}):
        self.gonio      = gonio
        if self.gonio.calc_psuedo == False:
            self.gonio.calc_psuedo = True
            self.gonio._update_psuedo()
        self.beam_slits = beam_slits
        self.det_slits  = det_slits
        self.sample     = sample
        # fraction horz polarization
        self.fh         = 1.0
        # arbitrary scale factor.
        # Note if Io = 1million cps
        # the using this scale makes the normalized intensity
        # close to cps
        self.scale      = 1.0e6

    ##########################################################################
    def ctot_stationary(self,plot=False):
        """
        correction factors for stationary measurements (e.g. images)
        """
        cp = self.polarization()
        cl = self.lorentz_stationary()
        ca = self.active_area(plot=plot)
        ct = (self.scale)*(cp)*(cl)*(ca)
        if plot == True:
            print "Polarization=%f" % cp
            print "Lorentz=%f" % cl
            print "Area=%f" % ca
            print "Scale=%f" % self.scale
            print "Total=%f" % ct
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
    def active_area(self,plot=False):
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
        alpha = self.gonio.pangles['alpha']
        beta  = self.gonio.pangles['beta']
        if plot == True:
            print 'alpha = ', alpha, ', beta = ', beta
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
                                     sample=sample,plot=plot)
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
    #print psic
    #
    beam_slits = {'horz':.6,'vert':.8}
    det_slits = {'horz':20.0,'vert':10.5}
    sample = {}
    sample['polygon'] = [[1.,1.], [.5,1.5], [-1.,1.], [-1.,-1.],[0.,.5],[1.,-1.]]
    sample['angles']  = {'phi':108.0007,'chi':0.4831}
    #
    cor = CtrCorrectionPsic(gonio=psic,beam_slits=beam_slits,
                            det_slits=det_slits,sample=sample)
    ct = cor.ctot_stationary(plot=True)

##########################################################################
if __name__ == "__main__":
    """
    test 
    """
    test1()
    #test2()


