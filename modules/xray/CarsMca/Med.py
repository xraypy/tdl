"""
Support for Multi-Element Detectors (Med).

Author:        Mark Rivers
Created:       Sept. 18, 2002.  Based on earlier IDL code.
Modifications:
"""

import Mca
import Numeric
import spline
import copy

#########################################################################
#class Med(Mca.Mca):
class Med:
    """
    The MED class is basically a collection of Mca objects.
    
    This class is device-independent.
    
    Its methods generally simply apply the Mca class methods to each Mca object
    in the collection. The Med class itself is most commonly used for reading
    data from disk files. More importantly, this class is the superclass of the
    epicsMed class.
    """

    def __init__(self, n_detectors=16, name='med'):
        """
        Initialization code for creating a new Med object.
    
        Keywords:
            n_detectors: The number of detectors (Mca objects) in the Med.

            name: Name of the mca object, eg a file name or detector name
        """
        self.n_detectors = n_detectors
        self.name        = name
        self.mcas        = []
        for i in range(n_detectors):
            mca_name = 'mca:%s' % str(i + 1)
            self.mcas.append(Mca.Mca(name=mca_name))
        self.environment = []
        # for deadtimes
        self.tau = None
        self.update = True

    def __repr__(self):
        lout = "Name = %s\n" % self.name
        lout = lout + "Number of detectors = %i\n" % self.n_detectors
        for mca in self.mcas:
            lout = lout + mca.__repr__()
        return lout

    def __copy__(self):
        new = Med()
        new.n_detectors = copy.copy(self.n_detectors)
        new.name        = copy.copy(self.name)
        new.mcas        = []
        for i in range(self.n_detectors):
            new.mcas.append(self.mcas[i].__copy__())
        return new
        
    def __deepcopy__(self,visit):
        new = Med()
        new.n_detectors = copy.copy(self.n_detectors)
        new.name        = copy.copy(self.name)
        new.mcas        = []
        for i in range(self.n_detectors):
            new.mcas.append(self.mcas[i].__deepcopy__(visit))
        return new

    ########################################################################
    def set_environment(self, environment):
        """
        Copies a list of McaEnvironment objects to the Med object.

        Inputs:
            environment:
                A list of McaEnvironment objects.
        """
        self.environment = environment

    #########################################################################
    def get_mcas(self):
        """
        Returns a list of Mca objects from the Med.
        """
        return self.mcas

    #########################################################################
    def initial_calibration(self, energy):
        """
        Performs an initial energy calibration for each Mca in the Med.

        Inputs:
            energy: The energy of the largest peak in the spectrum.
                
        See the documentation for Mca.initial_calibration() for more information.
        """
        for mca in self.mcas:
            mca.initial_calibration(energy)

    #########################################################################
    def final_calibration(self, peaks):
        """
        Performs a final energy calibration for each Mca in the Med.

        Inputs:
            peaks: A list of McaPeak objects. This list is typically read from a
                   disk file with function Mca.read_peaks().
                
        See the documentation for Mca.final_calibration() for more information.
        """
        print "Not yet implemented"
        #for mca in self.mcas:
        #    mca.final_calibration(peaks)

    #########################################################################
    def get_calibration(self):
        """
        Returns a list of McaCalibration objects, one for each Mca in the Med.
        """
        calibration = []
        for mca in self.mcas:
            calibration.append(mca.get_calibration())
        self.calibration = calibration
        return calibration

    #########################################################################
    def set_calibration(self, calibration):
        """
        This procedure sets the calibration parameters for the Med.
        The calibration information is contained in an object or list of 
        objects of type McaCalibration.
        
        Inputs:
            calibration: A single object or a list of objects of type McaCalibration
                         containing the calibration parameters for each Mca.
                         If a single object is passed then this is written to each Mca.
                         If a list of objects is passed then calibration[i] is written to
                         Mca[i].
        """
        if (isinstance(calibration, Mca.McaCalibration)):
            for mca in self.mcas:
                mca.set_calibration(calibration)
        else:  # Assume it is a list or tuple
            for i in range(self.n_detectors):
                self.mcas[i].set_calibration(calibration[i])


    #########################################################################
    def get_elapsed(self):
        """
        Returns the elapsed parameters for the Med.
        The elapsed information is contained in a list of structures of type
        McaElapsed.
        
        Outputs:
            Returns a list of structures of type McaElapsed.
            
        Procedure:
            This function simply invokes Mca.get_elapsed for each Mca in the Med
            and stores the results in the returned list.
        """
        elapsed = []
        for mca in self.mcas:
            elapsed.append(mca.get_elapsed())
        return elapsed

    #########################################################################
    def set_elapsed(self, elapsed):
        """
        Sets the elapsed parameters for the Med.
        The elapsed information is contained in an object or list of 
        objects of type McaElapsed.

        Inputs:
            elapsed: A single structure or a list of structures of type McaElapsed
                     containing the elapsed parameters for each Mca.
                     If a single object is passed then this is written to each Mca.
                     If a list of objects is passed then elapsed[i] is written to Mca[i].
        """
        if (isinstance(elapsed, Mca.McaElapsed)):
            for mca in self.mcas:
                mca.set_elapsed(elapsed)
        else:  # Assume it is a list or tuple
            for i in range(self.n_detectors):
                self.mcas[i].set_elapsed(elapsed[i])

    #########################################################################
    def get_rois(self,units="channel"):
        """
        Returns the region-of-interest information for each Mca in the Med.

        Keywords:
            units:  the returned ROI objects will be converted to the given units
                    valid values are "channel", "kev", "ev".  Default is "channel"

        Outputs:
            Returns a list of list of lists of McaRoi objects.
            The length of the outer list is self.n_detectors, the length of the
            list for each Mca is the number of ROIs defined for that Mca.
        """
        rois = []
        for mca in self.mcas:
            rois.append(mca.get_rois(units=units))
        return rois

    #########################################################################
    def get_roi_counts(self, background_width=1):
        """
        Returns the net and total counts for each Roi in each Mca in the Med.

        Outputs:
            Returns a tuple (total, net).  total and net are lists of lists
            containing the total and net counts in each ROI.  The length of the
            outer list is self.n_detectors, the length of the total and net lists
            list for each Mca is the number of ROIs defined for that Mca.
        """
        total = []
        net = []
        for mca in self.mcas:
            t, n = mca.get_roi_counts(background_width)
            total.append(t)
            net.append(n)
        return (total, net)

    #########################################################################
    def set_rois(self, rois):
        """
        This procedure sets the ROIs for the Med (blowing away old ones)
        The elapsed information is contained in a list of McaRoi objects,
        or list of such lists.

        Inputs:
            rois: A single list or a nested list of objects McaROI objects.
                  If a single list is passed then this is written to each Mca.
                  If a list of lists is passed then rois[i][*] is written to Mca[i].
        """
        if (len(rois) <= 1):
            for mca in self.mcas:
                mca.set_rois(rois)
        else:
            for i in range(self.n_detectors):
                self.mcas[i].set_rois(rois[i])

    #########################################################################
    def add_roi(self, roi):
        """
        This procedure adds an ROI to each Mca in the Med.

        Inputs:
            roi: A single McaROI to be added.
        """
        for mca in self.mcas:
            mca.add_roi(roi)
            
    #########################################################################
    def delete_roi(self, index):
        """
        This procedure deletes the ROI at position "index" from each Mca in the
        Med.

        Inputs:
            index:  The index number of the ROI to be deleted.
        """
        for mca in self.mcas:
            mca.delete_roi(index)

    #########################################################################
    def copy_rois(self, source_mca=0, energy=False):
        """
        This procedure copies the ROIs defined for one Mca in the Med to all of
        the other Mcas.

        Inputs:
            source_mca: The index number of the Mca from which the ROIs are to
                        be copied.  This number ranges from 0 to self.n_detectors-1.
                        The default is the first Mca (index=0).
                
        Keywords:
            energy: Set this keyword if the ROIs should be copied by their position
                    in energy rather than in channels. This is very useful when 
                    copying ROIs when the calibration parameters for each Mca in 
                    the Med are not identical.
        """
        if energy == True:
            units = "keV"
        else:
            units = "channel"
        rois = self.mcas[source_mca].get_rois(units=units)
        self.set_rois(rois)


    #########################################################################
    def update_correction(self,tau=None):
        """ Update mca deadtime correction factors """
        if tau is not None:
            if len(tau) != self.n_detectors:
                print "Error: tau array must be of length %d" % self.n_detectors
                self.tau = None
            else:
                self.tau = tau
                self.update = True
        if self.update == True:
            for j in range(self.n_detectors):
                if self.tau != None:
                    t = self.tau[j]
                else:
                    t = None
                self.mcas[j].update_correction(tau=t)
            self.update = False

    #########################################################################
    def get_data(self, detectors=[], total=False, align=False, correct=True):
        """
        Returns the data from each Mca in the Med as a 2-D Numeric array
        
        Keywords:
            total: Set this keyword to return the sum of the spectra from all
                   of the Mcas as a 1-D Numeric array.
                
            align: Set this keyword to return spectra which have been shifted and
                   and stretched to match the energy calibration parameters of the
                   first detector.  This permits doing arithmetic on a
                   "channel-by-channel" basis. This keyword can be used alone
                   or together with the TOTAL keyword, in which case the data
                   are aligned before summing.

            detectors: A list of detectors to be passed back.  An empty
                       list (default) means use all the detectors.
                       Note detector indexing starts at zero!

            correct:
                True means to apply deadtime correction, false ignores it
                
        Outputs:
            By default this function returns a long 2-D array of counts dimensioned
            [self.n_detectors,nchans]
            
            If the "total" keyword is set then the array dimensions are [1,nchans]
        """
        
        if len(detectors) > 0 and len(detectors) < self.n_detectors:
            det_idx = detectors
        else:
            det_idx = range(self.n_detectors)
            
        temp = self.mcas[det_idx[0]].get_data()
        nchans = len(temp)

        #data = Numeric.zeros((self.n_detectors, nchans))
        data = Numeric.zeros((len(det_idx), nchans))

        #for i in range(self.n_detectors):
        for i in range(len(det_idx)):
            d = int(det_idx[i])
            data[i,:] = self.mcas[d].get_data(correct=correct)

        if align ==True and len(det_idx) > 1:
            #ref_energy = self.mcas[0].get_energy()
            d = int(det_idx[0])
            ref_energy = self.mcas[d].get_energy()
            for i in range(len(det_idx)):
                d = int(det_idx[i])
                energy = self.mcas[d].get_energy()
                temp = spline.spline_interpolate(energy, data[i,:], ref_energy)
                data[i,:] = (temp+.5).astype(Numeric.Int)

        if total == True and len(det_idx) > 1:
            d = Numeric.sum(data)
            return [d]
        else:
            return data

    #########################################################################
    def get_energy(self,detectors = []):
        """
        Returns a list of energy arrays, one array for each Mca in the Med.
        See the documentation for Mca.get_energy() for more information.
        """
        if len(detectors) > 0 and len(detectors) < self.n_detectors:
            det_idx = detectors
        else:
            det_idx = range(self.n_detectors)

        energy = []
        #for mca in self.mcas:
        for i in range(len(det_idx)):
            energy.append(self.mcas[int(det_idx[i])].get_energy())
        return energy
