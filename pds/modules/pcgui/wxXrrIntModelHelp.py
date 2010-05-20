########################################################################
"""
Tom Trainor (fftpt@uaf.edu)
This should be run as a child window from wxXrrBuilder

--------------
Modifications:
--------------

"""
########################################################################

from PythonCard import model, dialog
import wx
import os
import types

from   wxUtil import wxUtil

########################################################################

XrrIntModelParamHelp = """ XRR Interface Model Help

This window provides access to the xrr.interface_model module,
used for computing various interfacial chemical distributions and
the resulting reflectivity and fluorescent yield.
To construct a model see the XRR model builder window. 

==== Group/Model ====
** Use the update button to refresh the list of pds variables...

** Select an InterfaceModel instance as the input model
   (see Model Builder to initiate an instance)

==== Components / Component Distribution Options ====
** Select a component defined within the model to specify its interface
   distribution.
   
   Note that the component concentration defined in the top and bottom
   layers of the model is fixed and cannot be varied (see model builder
   to change these values).  These values are displayed in the boxes labeled
   'Top' and 'Subs'.

** The Norm flag has the following allowed values
   (and applies only to the selected component):
   0 = No normalization
   1 = Normalize to the original model interface mass balance 
   2 = Normalize so that CX[1]    = CX[0]    (CX[0] = substrate concentration)
   3 = Normalize so that CX[nz-2] = CX[nz-2] (CX[nz-2] = top concentration)

** If the components 'Scale to Norm' check box is selected, the model will
   apply the normalization factor (that satifies the above condition)
   to all the distribution scale (CX) parameters for the selected component.  

** The Density Norm flag has the following allowed values
   (and applies globally)
   0 = No normalization
   1 = Normalize such that rho[1]    = rho[0]    (rho[0] = substrate density)
   2 = Normalize such that rho[nz-2] = rho[nz-1] (rho[nz-1] = top density)

** If the Density 'Scale to Norm' check box is selected, the model will
   apply the normalization factor (that satifies the above condition)
   to all the distribution scale (CX) parameters.  This is global (ie
   applies to all components).

   Note when using density normalization it is applied after any calculating
   the component distribution (and applying any component normalization).
   Therefore, the density scale factor necessary is computed from the
   initially calculated distribution and then applied to recompute the scaled
   component/element distributions.  

** If the 'Auto Calc Dist' button is selected, the distribution function
   will be recomputed whenever a parameter changes.  Otherwise you must
   hit the 'Calc Dist' button

** Multiple distribution models may be defined in the interface region.
   The final component distribution is the sum over all defined distributions.
   (note normalization is applied to the final sum)

** The interface region always starts at Z = 0 and ends at the bottom
   of the top layer (therefore valid Z-ranges are from 0 to zst_top)

** For most distribution models using Zstart and Zend values of 'None' means
   allow the model to extend throughout the entire interface region.

==== Distribution Models ====

** Box model (box):
   - CX = concentration in moles/cm^3
   - This fixed concentration extends from Zstart to Zend
     (if Zstart = None, then model start at Z=0)
     (if Zend = None, then model ends at top of the interface)
     
** Linear model (linear):
   - CX = left hand side concetration in moles/cm^3
   - CXen = right hand concentration in moles/cm^3
   - The concentration varies linearly from Zstart to Zend 
     (if Zstart = None, then model start at Z=0)
     (if Zend = None, then model ends at top of the interface)

** Error Function (erf/erfc):
   - CX = max concentration in moles/cm^3
   - Center = Z-value for 1/2 CX (angstroms)
   - Width = Width of the roll over (angstroms).
   - This model can effectively reproduce a half box model using
     Width = 0.
   - The erfc model is the complement, ie the plateau is on the
     high Z side.  

** Exponential Function (exp/expc):
   - CX = max concentration in moles/cm^3
   - Center = Z-value where transition from CX to exp form (angstroms)
   - Width = 1/e length (angstroms).
   - The expc model is the complement, ie the plateau is on the
     high Z side.  

** Gaussian (gauss):
   - Zstart and Zend are ignored (calculation is done for entire interface)
   - CX = Gaussian amplitude in moles/cm^3. 
   - Center = Gaussian center (angstroms)
   - Width = Gaussian standard deviation (angstroms).
   - Note the function is normalized so that the max = CX

==== Data ====
** Select arrays for theta, Rdata and FYdata
     
==== Params ====

Note to change a parameter in the model, adjust its value then click
update paramters.  i.e. unlike the distribution parameters these are
not automatically updated.  

** Energy = Incident beam energy (eV)

** Convolution width = width of convolution function (in degrees)

** Sample length = length of the sample (mm)

** Beam vert =  vert beam size (mm)

** Beam horz = horz beam size (mm)

** Model Area variation = Area flag to indicate if area calculation
   should be performed.  ie model predicts FY change due to acive area
   variation with theta: 0.0 no, 1.0 yes

** FY Element = Element symbol or Z value for FY calcs
   Use 0 to turn off the FY calculations

** FY Energy =  Energy (eV) of FY line (e.g. 'Fe Ka')

** FY Detector Angle = detector angle, ie take off angle btwn
   substrate and det. (deg)

** Theta Norm =  theta value to use for yield normalization. (deg)

** Roughness Flag = roughness flag
   0.0 ignore roughness effect on t,
   1.0 use t = 1+r*exp(-(q*sigma)^2)
   Note for slabified model using slab delta > 0 (see model builder)
   the roughness values within the interface region are ignored!  ie
   we assume that roughness is modeled with density variations...

** Integration Delta Z = delta z for FY int (angstroms)

** FY Base base penetration depth factor = Mulitple of the
   substrate penetration depth to use as integration
   depth.  Use 0.0 to ignore, ie use d[0] (given substrate thickness)
"""

###########################################################
class wxXrrIntModelHelp(model.Background, wxUtil):

    ###########################################################
    def on_initialize(self, event):
        # Initialization
        # including sizer setup, do it here
        # self.setupSizers()
        self.help = XrrIntModelParamHelp
        self.components.HelpText.text = self.help

    ###########################################################
    def on_menuFileExit_select(self,event): 
        self.close()
