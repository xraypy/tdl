"""
Wrapper for the reflectivity / reflection XSW functions in _ref.dll

Authors/Modifications:
----------------------
T. Trainor (tptrainor@alaska.edu)

Notes:
------
* For this module to work the directory 'wrap/libs'
  must be in your shell's execution path so the load_library
  command can find _ref.dll (which also needs to find the gsl dll's)

"""
#######################################################################

import ctypes as C
import numpy as num
import os

#######################################################################

#Get the name of the 'libs' directory under wrap
libspath = os.path.dirname(__file__)
libspath = os.path.join(libspath,'libs')
libspath = os.path.abspath(libspath)

# import the dll 
_ref_dll = num.ctypeslib.load_library('_ref',libspath)

#######################################################################
_ref_dll._Reflectivity.restype = C.c_int

argtypes = [C.c_int,                  # nlayer
            C.c_int,                  # nelem
            C.c_int,                  # nthet
            C.POINTER(C.c_double),    # calc_params
            #
            C.POINTER(C.c_double),    # d
            C.POINTER(C.c_double),    # rho
            C.POINTER(C.c_double),    # sigma
            C.POINTER(C.POINTER(C.c_double)), # comp ??
            #
            C.POINTER(C.c_double),    # elem_z
            C.POINTER(C.c_double),    # fp
            C.POINTER(C.c_double),    # fpp
            C.POINTER(C.c_double),    # amu
            C.POINTER(C.c_double),    # mu_at
            #
            C.POINTER(C.c_double),    # theta_arr
            C.POINTER(C.c_double),    # R
            C.POINTER(C.c_double),    # Y
            #
            C.POINTER(C.c_double),    # delta
            C.POINTER(C.c_double),    # bet
            C.POINTER(C.c_double),    # amu_t
            C.POINTER(C.c_double),    # mu_t
            #
            C.POINTER(C.c_double),    # Re_X
            C.POINTER(C.c_double),    # Im_X
            C.POINTER(C.c_double),    # Re_Ai
            C.POINTER(C.c_double),    # Im_Ai
            C.POINTER(C.c_double),    # Re_Ar
            C.POINTER(C.c_double),    # Im_Ar
            C.POINTER(C.c_double),    # Re_g
            C.POINTER(C.c_double)]    # Im_g
            
# not sure how to deal with comp array in argtypes? Its a double **
#_ref_dll._Reflectivity.argtypes = argtypes

#######################################################################
class _Reflectivity:
    """
    simple class to interface with c-library reflectivity calc
    """
    def _init_model(self, nlayer, nelem, nthet):
        """
        create 'empty' model arrays. This method can/should be
        bypassed by other model creation routines....
        """
        #
        self.nlayer      = int(nlayer)
        self.nelem       = int(nelem)
        self.nthet       = int(nthet)
        #
        self._init_params()
        #
        self.d           = num.zeros(self.nlayer, dtype=num.double)
        self.rho         = num.zeros(self.nlayer, dtype=num.double)
        self.sigma       = num.zeros(self.nlayer-1, dtype=num.double)
        self.comp        = num.zeros((self.nelem, self.nlayer), dtype=num.double )
        #
        self.elem_z      = num.zeros(self.nelem, dtype=num.double)
        self.fp          = num.zeros(self.nelem, dtype=num.double)
        self.fpp         = num.zeros(self.nelem, dtype=num.double)
        self.amu         = num.zeros(self.nelem, dtype=num.double)
        self.mu_at       = num.zeros(self.nelem, dtype=num.double)
        #
        self.theta       = num.zeros(self.nthet, dtype=num.double)

    ###################################################################
    def _init_params(self,):
        """
        # calc_params[0] = energy (eV)
        # calc_params[1] = width of convolution function (in degrees)
        # calc_params[2] = length of the sample (mm)
        # calc_params[3] = vert beam size (mm)
        # calc_params[4] = horz beam size (mm)
        # calc_params[5] = Area flag to indicate if area calculation should be performed
        #                  (ie model predicts FY change due to acive area variation with theta)
        #                  0.0 no, 1.0 yes
        # calc_params[6] = fy_idx
        # calc_params[7] = fy energy (eV)
        # calc_params[8] = detector angle (ie take off angle btwn substrate and det) (deg)
        # calc_params[9] = theta norm (theta value to use for yield normalization) (deg)
        # calc_params[10] = roughness flag
        #                  0.0 ignore roughness effect on t,
        #                  1.0 use t = 1+r*exp(-(q*sigma)^2)
        # calc_params[11] = delta z for FY int (angstroms)
        # calc_params[12] = base penetration depth factor
        #                   0.0 (or less than zero) means ignore, ie use d[0]
        # calc_params[13] = reflectivity scale factor
        #
        """
        self.calc_params = num.zeros(14, dtype=num.double)
        self.calc_params[0] = 10000.  # energy (eV)
        self.calc_params[1] = 0.01    # wconv (deg)
        self.calc_params[2] = 50.0    # sample len (mm)
        self.calc_params[3] = 0.05    # bvert (mm)
        self.calc_params[4] = 10.0    # bhorz (mm)
        self.calc_params[5] = 0.0     # aflag
        self.calc_params[6] = 0.0     # fy_idx
        self.calc_params[7] = 0.0     # fy energy (eV)
        self.calc_params[8] = 90.0    # detang (deg) 
        self.calc_params[9] = 1.0     # tnorm (deg)
        self.calc_params[10] = 1.0    # rflag
        self.calc_params[11] = 10.0   # del z (ang)
        self.calc_params[12] = 3.0    # pdeth
        self.calc_params[13] = 1.0    # rscale
        

    ###################################################################
    def _init_calc_arrays(self):
        """
        create/re-init calc arrays. These arrays are filled in
        by the _Reflectivity function in C.  They must be correctly
        sized, but can be inited to zeros...
        """
        #
        self.R           = num.zeros(self.nthet, dtype=num.double)
        self.Y           = num.zeros(self.nthet, dtype=num.double)
        #
        self.delta       = num.zeros(self.nlayer, dtype=num.double)
        self.beta        = num.zeros(self.nlayer, dtype=num.double)
        self.amu_t       = num.zeros(self.nlayer, dtype=num.double)
        self.mu_t        = num.zeros(self.nlayer, dtype=num.double)
        #
        self.Re_X        = num.zeros(self.nlayer, dtype=num.double)
        self.Im_X        = num.zeros(self.nlayer, dtype=num.double)
        self.Re_Ai       = num.zeros(self.nlayer, dtype=num.double)
        self.Im_Ai       = num.zeros(self.nlayer, dtype=num.double)
        self.Re_Ar       = num.zeros(self.nlayer, dtype=num.double)
        self.Im_Ar       = num.zeros(self.nlayer, dtype=num.double)
        self.Re_g        = num.zeros(self.nlayer, dtype=num.double)
        self.Im_g        = num.zeros(self.nlayer, dtype=num.double)

    ###################################################################
    def _arr_ptrs(self):
        """
        get pointers to pass to c-function
        """
        dptr = C.POINTER(C.c_double)
        #
        self.calc_params_ptr = self.calc_params.ctypes.data_as(dptr)
        #
        self.d_ptr           = self.d.ctypes.data_as(dptr)
        self.rho_ptr         = self.rho.ctypes.data_as(dptr)
        self.sigma_ptr       = self.sigma.ctypes.data_as(dptr)
        self.comp_ptr        = (dptr*self.nelem)(*[row.ctypes.data_as(dptr)
                                                   for row in self.comp])
        #
        self.elem_z_ptr      = self.elem_z.ctypes.data_as(dptr)
        self.fp_ptr          = self.fp.ctypes.data_as(dptr)
        self.fpp_ptr         = self.fpp.ctypes.data_as(dptr)
        self.amu_ptr         = self.amu.ctypes.data_as(dptr)
        self.mu_at_ptr       = self.mu_at.ctypes.data_as(dptr)
        #
        self.theta_ptr       = self.theta.ctypes.data_as(dptr)
        self.R_ptr           = self.R.ctypes.data_as(dptr)
        self.Y_ptr           = self.Y.ctypes.data_as(dptr)
        #
        self.delta_ptr       = self.delta.ctypes.data_as(dptr)
        self.beta_ptr        = self.beta.ctypes.data_as(dptr)
        self.amu_t_ptr       = self.amu_t.ctypes.data_as(dptr)
        self.mu_t_ptr        = self.mu_t.ctypes.data_as(dptr)
        #
        self.Re_X_ptr        = self.Re_X.ctypes.data_as(dptr)
        self.Im_X_ptr        = self.Im_X.ctypes.data_as(dptr)
        self.Re_Ai_ptr       = self.Re_Ai.ctypes.data_as(dptr)
        self.Im_Ai_ptr       = self.Im_Ai.ctypes.data_as(dptr)
        self.Re_Ar_ptr       = self.Re_Ar.ctypes.data_as(dptr)
        self.Im_Ar_ptr       = self.Im_Ar.ctypes.data_as(dptr)
        self.Re_g_ptr        = self.Re_g.ctypes.data_as(dptr)
        self.Im_g_ptr        = self.Im_g.ctypes.data_as(dptr) 

    ###################################################################
    def _calc(self,init_ptrs=True,init_arrs=False):
        """
        Calculate reflectivity/yield
        The C function has the following call:
        _Reflectivity(
            int nlayer, int nelem, int nthet, double *calc_params,  
            double *d, double *rho, double *sigma, double **comp,
            double *elem_z, double *fp, double *fpp, double *amu, double *mu_at,
            double *theta_arr, double *R, double *Y, 
            double *del, double *bet, double *amu_t, double *mu_t,
            double *Re_X, double *Im_X, double *Re_Ai, double *Im_Ai, 
            double *Re_Ar, double *Im_Ar, double *Re_g, double *Im_g )
        """
        if init_arrs: self._init_calc_arrays()
        if init_ptrs: self._arr_ptrs()
        
        ret = _ref_dll._Reflectivity( C.c_int(self.nlayer),
                                      C.c_int(self.nelem),
                                      C.c_int(self.nthet),
                                      self.calc_params_ptr,
                                      #
                                      self.d_ptr,
                                      self.rho_ptr,
                                      self.sigma_ptr,
                                      self.comp_ptr,
                                      #
                                      self.elem_z_ptr,
                                      self.fp_ptr,
                                      self.fpp_ptr,
                                      self.amu_ptr,
                                      self.mu_at_ptr,
                                      #
                                      self.theta_ptr,
                                      self.R_ptr,
                                      self.Y_ptr,
                                      #
                                      self.delta_ptr,
                                      self.beta_ptr,
                                      self.amu_t_ptr,
                                      self.mu_t_ptr,
                                      #
                                      self.Re_X_ptr,
                                      self.Im_X_ptr,
                                      self.Re_Ai_ptr,
                                      self.Im_Ai_ptr,
                                      self.Re_Ar_ptr,
                                      self.Im_Ar_ptr,
                                      self.Re_g_ptr,
                                      self.Im_g_ptr )

    ###################################################################
    def _show_ptrs(self):
        print 'nlayer      ', C.c_int(self.nlayer)
        print 'nelem       ', C.c_int(self.nelem)
        print 'nthet       ', C.c_int(self.nthet)
        print 'calc params ', self.calc_params_ptr
        #
        print 'd           ', self.d_ptr
        print 'rho         ', self.rho_ptr
        print 'sigma       ', self.sigma_ptr
        print 'comp        ', self.comp_ptr
        #
        print 'elem z      ', self.elem_z_ptr
        print 'fp          ', self.fp_ptr
        print 'fpp         ', self.fpp_ptr
        print 'amu         ', self.amu_ptr
        print 'mu_at       ', self.mu_at_ptr
        #
        print 'theta       ', self.theta_ptr
        print 'R           ', self.R_ptr
        print 'Y           ', self.Y_ptr
        #
        print 'delta       ', self.delta_ptr
        print 'beta        ', self.beta_ptr
        print 'amu_t       ', self.amu_t_ptr
        print 'mu_t        ', self.mu_t_ptr
        #
        print 'Re_X        ', self.Re_X_ptr
        print 'Im_X        ', self.Im_X_ptr
        print 'Re_Ai       ', self.Re_Ai_ptr
        print 'Im_Ai       ', self.Im_Ai_ptr
        print 'Re_Ar       ', self.Re_Ar_ptr
        print 'Im_Ar       ', self.Im_Ar_ptr
        print 'Re_g        ', self.Re_g_ptr
        print 'Im_g        ', self.Im_g_ptr

#######################################################################
#######################################################################
def lib_test():
    x = num.array([[1.1,2.1,3.1],[4.1,5.1,6.1]])
    y = num.array([-1.,-20.,-30.])
    dptr = C.POINTER(C.c_double)
    xptr = (dptr*len(x))(*[row.ctypes.data_as(dptr) for row in x])
    yptr = y.ctypes.data_as(dptr)
    nr   = len(x)
    nc   = len(x[0])
    _ref_dll._test(xptr,yptr,nr,nc)

#######################################################################
if __name__ == "__main__":
    #lib_test()
    print dir()
    
