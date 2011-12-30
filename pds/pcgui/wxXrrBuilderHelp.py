"""
This should be run as a child window from wxXrrBuilder

"""
########################################################################

from PythonCard import model, dialog
import wx
import os
import types

from  pds.pcgui.wxUtil import wxUtil

########################################################################

XrrBuilderParamHelp = """ XRR Help

This window provides an interface to construct reflectivity models
and computing reflectivity and fluorescent yield profiles.
Various interfacial chemical distributions may be computed using the
interface_model applications.

==== Model ====
** Read in an existing model using the 'Model' menu

** Provide a model name to save/write to in the Save Model field

** Slab Delta Z = target slab thickness for 'slabified' model
   If this value is less than or equal to zero, the layer thicknesses
   are used as the to construct the slab model (ie slab model = layer model)
   
   Note for slabified model using slab Delta > 0 the roughness values
   within the interface region are ignored!  ie we assume that roughness
   is modeled with density variations (see InterfaceModel module)

==== Layers ====
** Composition = layer composition in terms of chemical components.
   Use the following syntax (eg):
     {(1.0)[Si_O2]}{(1e-7)[Fe_O]}{(0.12)[Ca0.8_Na0.4_O1]}

** Thickness = layer thickness in angstroms.

** Density = layer density in g/cm^3.

** Roughness = layer roughness in angstroms.

==== Data ====
** Select arrays for theta, Rdata and FYdata

==== Params ====

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
class wxXrrBuilderHelp(model.Background, wxUtil):

    ###########################################################
    def on_initialize(self, event):
        # Initialization
        # including sizer setup, do it here
        # self.setupSizers()
        self.help = XrrBuilderParamHelp
        self.components.HelpText.text = self.help

    ###########################################################
    def on_menuFileExit_select(self,event): 
        self.close()
