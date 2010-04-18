"""
Function for extracting RASD data from spec "raxr" scans

Authors/Modifications:
----------------------
Frank Heberling (Frank.Heberling@kit.edu)

Notes:
------
Currently probably specific for APS 13ID-C

To use the function you may type in these 4 variables with the values you need:

    spec_path = 'D:/Data/Surf_Diff/Aps_Data_ID' (path to your spec file)
    spec = 'caco3_july09h.spc'                  (name of your spec file)
    first_scan = 45 
    last_scan = 143 (first and last hklscan in your spec file, that belong to the raxr scan that is to be integrated)
    
and run it like this:

    A = scandata.rasd.rasd(spec_path, spec, first_scan, last_scan)

"rasd" appends all the hklscans specified to one RASDData object (the two images and the 'io's in one
hklscan are summed up).
After interactively integrating the images in this object (as a ScanData object) you'll end up with a plot
of your RASD-wiggle, a .rsd file (named rasdHK_L.rsd with HKL being the values of the corresponding scan)
that contains the plotted data, and a RASDData object called "A" that contains:

A.H, A.K, A.L (single values from first hklscan)
A.E (array of energy values)
A.image (list of images)
A.image_peaks (dictionary with extracted intensity and background values) 
A.scalers (dictionary with exact H,K,L values, io, and beta values)
A.positioners (dictionary with diffractometer angles)
A.G (state['G'] matrix from first hklscan)
A.corr (array of values: corr = Ci/(Cp * I0))
A.F, A.Ferr (arrays with structure factors / uncertainty of structure factors, respectively)

Todo:
-----
- make RASDData use the more sophisticated correction functions from CtrData
- include active area correction
- include normalization, amplitude, and phase fitting
- Test

"""
#############################################################################

import types
import pylab
import numpy as Num
from   reader import Reader
from   data import ScanData
import image_data
from image_menu import image_menu

#############################################################################
def rasd(spec_path,spec,first_scan,last_scan):
    """
    Merges all the hklscans of a spec "raxr" scan to one RASDData object containing images and info from the
    spec file, integrates the images interactively (as a ScanData object), calculates polarization + lorentz correction
    factors, plots the RASD wiggle, and writes the Information into a .rsd file.
    """
    spec_file = Reader(spec, spec_path)
    spec_file.image = False
    nscans = last_scan - first_scan + 1

    S = []
    for i in range(nscans):
        tmp = spec_file.spec_scan(first_scan + i, image = True)
        S.append(tmp)

    Z = RASDData()
    Z.append_rasd(S)
    data = ScanData('',[len(Z.E)],Z.scalers,Z.positioners,primary_axis=[],primary_det=None,state={'G':Z.G},
                    med=None,xrf=None,xrf_lines=None,image=Z.image,image_rois=None)
    image_menu(data)
    Z.image_peaks = data.image.peaks
    Z.corr = Z.calc_correction(data)
    Z.F    = (Z.corr*Z.image_peaks['I'])**0.5
    Z.Ferr = (Z.corr*Z.image_peaks['Ierr'])**0.5
    Z.dims = len(Z.E)

    Z.plot()
    filename = 'rasd'+str(int(Z.H))+str(int(Z.K))+'_'+str(Z.L)+'.rsd'
    Z.write_rsd(filename)
    return Z

#############################################################################
class RasdData:
    def __init__(self,scans=[]):
        #
        self.H    = 0
        self.K    = 0
        self.L    = 0
        self.dims = 0
        self.G    = Num.array([],float)
        #
        self.E     = Num.array([],float)
        self.image = []
        self.image_peaks = {}
        #
        self.scalers = {'H': Num.array([]),
                        'K': Num.array([]),
                        'L': Num.array([]),
                        'io': Num.array([]),
                        'Beta': Num.array([])}
        
        self.positioners = {'nu': Num.array([]),
                            'eta': Num.array([]),
                            'del': Num.array([]),
                            'mu': Num.array([]),
                            'phi': Num.array([]),
                            'chi': Num.array([])}

        self.corr  = Num.array([],float)
        self.F     = Num.array([],float)
        self.Ferr  = Num.array([],float)

    ##############################################################################
    def append_rasd(self, scans):
        self.H    = scans[0].scalers['H'][0]
        self.K    = scans[0].scalers['K'][0]
        self.L    = scans[0].scalers['L'][0]
        self.G    = scans[0].state['G']
        for scan in scans:
            self.E = Num.append(self.E, scan['ENERGY'][0])
            self.image.append(scan.image.image[0] + scan.image.image[1])
            self.scalers['H']    = Num.append(self.scalers['H'],scan.scalers['H'][0])
            self.scalers['K']    = Num.append(self.scalers['K'],scan.scalers['K'][0])
            self.scalers['L']    = Num.append(self.scalers['L'],scan.scalers['L'][0])
            self.scalers['io']   = Num.append(self.scalers['io'],scan.scalers['io'][0]+scan.scalers['io'][1])
            self.scalers['Beta'] = Num.append(self.scalers['Beta'],scan.scalers['Beta'][0])
            self.positioners['phi']  = Num.append(self.positioners['phi'],scan.positioners['phi'][0])
            self.positioners['chi']  = Num.append(self.positioners['chi'],scan.positioners['chi'][0])
            self.positioners['mu']   = Num.append(self.positioners['mu'],scan.positioners['mu'][0])
            self.positioners['nu']   = Num.append(self.positioners['nu'],scan.positioners['nu'][0])
            self.positioners['eta']  = Num.append(self.positioners['eta'],scan.positioners['eta'][0])
            self.positioners['del']  = Num.append(self.positioners['del'],scan.positioners['del'][0])            

    ##############################################################################            
    def plot(self):
        """
        plot the rasd-wiggle
        """
        hkl = "%i   %i   %5.3f" % (self.H,self.K,self.L)
        figtitle = 'RASD scan at (hkl): '+hkl
        pylab.figure()
        pylab.plot(self.E, self.F,'b')
        pylab.title(figtitle)
        pylab.errorbar(self.E,self.F,self.Ferr, fmt ='o')
        pylab.xlabel('Energy (eV)')
        pylab.ylabel('F (arb. units)')
        pylab.show()

    ##########################################################################
    def write_rsd(self,fname = 'rasd.rsd'):
        """
        write data file
        """
        f = open(fname, 'w')
        f.write('#idx     Energy            F           F_err \n')
        for i in range(len(self.E)):
            if self.F[i] > 0:
                line = "%i   %6.3f    %6.6g    %6.6g\n" % (i,
                                                           self.E[i],
                                                           self.F[i],
                                                           self.Ferr[i])
                f.write(line)
        f.close()

    ###########################################################################       
    def calc_correction(self, data):
        """
        correction factors for pilatus
        """
        #calculate correction factors for specular rods
        if (Num.abs(self.H) < 0.01 and \
            Num.abs(self.K) < 0.01) :
            Cp = 1 - Num.sin(Num.radians(self.positioners['nu']))**2
            Ci = Num.sin(Num.radians(self.positioners['nu']/2))
        #or calculate correction factors for off-specular rods
        else:
            Cp = 1 - Num.cos(Num.radians(self.positioners['del']))**2 * \
                 Num.sin(Num.radians(self.positioners['nu']))**2
            Ci = Num.cos(Num.radians(self.positioners['del'])) * \
                 Num.sin(Num.radians(self.scalers['Beta']))

        corr = Ci/(Cp* self.scalers['io'])
        return corr
