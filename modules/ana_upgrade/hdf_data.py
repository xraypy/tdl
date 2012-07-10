"""
Scattering data

Authors / Modifications:
------------------------
Craig and Tom

Notes:
------

      
Todo:
-----
* HdfDataFile also needs search / filter methods
that return lists of point numbers


"""
##############################################################################
import numpy
import h5py

from tdl.modules.ana_upgrade import file_locker
from tdl.modules.ana_upgrade import image_data

##############################################################################

# Standard data attributes for a point (ie read returns, write requires)
DEFAULT_DETECTOR = {'name':'default',
                    'detector_params':{},
                    'integrate_params':{},
                    'correction_params':{},
                    'results':{}}
"""
note for image detector detector_params are:
    (yaw, pitch, roll, dist, cen_pixel, del_x, del_y, scale)
for correction:
    (bad_pixels, flat_files, spatial)    
"""
DEFAULT_POINT_DATA = {'name':'default',
                      'type':'',
                      'info':'',
                      'geom':'',
                      'hist':'',
                      'date':0,
                      'geom_angle_lbls':[],
                      'geom_angles':[],
                      'energy':0,
                      'lattice':[],
                      'q':0,
                      'orient_lbls':[],
                      'orient':[],
                      'vref_lbls':[],
                      'vref':[],
                      'position_lbls':[],
                      'positions':[],
                      'scaler_lbls':[],
                      'scaler':[],
                      'scaler_scale':[],
                      'detector_1':DEFAULT_DETECTOR}

# DOES NOT INCLUDE THE VALUES IN 'position_values'
# OR 'scaler_labels' BECAUSE THESE CHANGE BASED ON
# THE SCAN TYPE
GEN_KEYS = {'chi': ['angle_values', 0],
            'del': ['angle_values', 1],
            'eta': ['angle_values', 2],
            'mu': ['angle_values', 3],
            'nu': ['angle_values', 4],
            'phi': ['angle_values', 5],
            'real_a': ['lattice_values', 0],
            'real_b': ['lattice_values', 1],
            'real_c': ['lattice_values', 2],
            'real_alpha': ['lattice_values', 3],
            'real_beta': ['lattice_values', 4],
            'real_gamma': ['lattice_values', 5],
            'recip_a': ['lattice_values', 6],
            'recip_b': ['lattice_values', 7],
            'recip_c': ['lattice_values', 8],
            'recip_alpha': ['lattice_values', 9],
            'recip_beta': ['lattice_values', 10],
            'recip_gamma': ['lattice_values', 11],
            'lambda': ['lattice_values', 12],
            'or0_h': ['or_values', 0],
            'or0_k': ['or_values', 1],
            'or0_L': ['or_values', 2],
            'or0_del': ['or_values', 3],
            'or0_eta': ['or_values', 4],
            'or0_chi': ['or_values', 5],
            'or0_phi': ['or_values', 6],
            'or0_nu': ['or_values', 7],
            'or0_mu': ['or_values', 8],
            'or0_lambda': ['or_values', 9],
            'or1_h': ['or_values', 10],
            'or1_k': ['or_values', 11],
            'or1_L': ['or_values', 12],
            'or1_del': ['or_values', 13],
            'or1_eta': ['or_values', 14],
            'or1_chi': ['or_values', 15],
            'or1_phi': ['or_values', 16],
            'or1_nu': ['or_values', 17],
            'or1_mu': ['or_values', 18],
            'or1_lambda': ['or_values', 19]}#,
            #'Psi': ['position_values', 0],
            #'H': ['position_values', 1],
            #'K': ['position_values', 2],
            #'L': ['position_values', 3],
            #'Alpha': ['position_values', 4],
            #'Beta': ['position_values', 5],
            #'Epoch': ['scaler_values', 0],
            #'Seconds': ['scaler_values', 1],
            #'i1': ['scaler_values', 2],
            #'Bicron': ['scaler_values', 3],
            #'AmpTek_sc': ['scaler_values', 4],
            #'ROI1': ['scaler_values', 5],
            #'ROI2': ['scaler_values', 6],
            #'ROI3': ['scaler_values', 7],
            #'io': ['scaler_values', 8],
            #'IROI': ['scaler_values', 9]}

ATT_KEYS = ['date_stamp',
            'energy',
            'geom',
            'hist.%i',
            'info',
            'name',
            'type']

MISC_KEYS = ['Q',
             'haz']

# DO NOT INCLUDE 'image_data' OR 'corrected_image'
# AS THIS WILL CAUSE THEM TO BE (OVER)WRITTEN IN
# THE FILE.
DET_KEYS = {'bad_pixel_map': ['det_%i/corr_values.%i', 0],
            'bad_point': ['det_%i/corr_values.%i', 1],
            'image_changed': ['det_%i/corr_values.%i', 2],
            'image_max': ['det_%i/corr_values.%i', 3],
            'pixel_map_changed': ['det_%i/corr_values.%i', 4],
            'real_image_max': ['det_%i/corr_values.%i', 5],
            'rotangle': ['det_%i/corr_values.%i', 6],
            'sample_angles': ['det_%i/corr_values.%i', 7],
            'sample_diameter': ['det_%i/corr_values.%i', 8],
            'sample_polygon': ['det_%i/corr_values.%i', 9],
            'scale': ['det_%i/corr_values.%i', 10],
            'beam_slits': ['det_%i/det_values.%i', 0],
            'det_slits': ['det_%i/det_values.%i', 1],
            'bgrflag': ['det_%i/int_values.%i', 0],
            'cnbgr': ['det_%i/int_values.%i', 1],
            'compress': ['det_%i/int_values.%i', 2],
            'cpow': ['det_%i/int_values.%i', 3],
            'ctan': ['det_%i/int_values.%i', 4],
            'cwidth': ['det_%i/int_values.%i', 5],
            'filter': ['det_%i/int_values.%i', 6],
            'integrated': ['det_%i/int_values.%i', 7],
            'nline': ['det_%i/int_values.%i', 8],
            'rnbgr': ['det_%i/int_values.%i', 9],
            'roi': ['det_%i/int_values.%i', 10],
            'rpow': ['det_%i/int_values.%i', 11],
            'rtan': ['det_%i/int_values.%i', 12],
            'rwidth': ['det_%i/int_values.%i', 13],
            'alpha': ['det_%i/result_values.%i', 0],
            'beta': ['det_%i/result_values.%i', 1],
            'ctot': ['det_%i/result_values.%i', 2],
            'F': ['det_%i/result_values.%i', 3],
            'F_changed': ['det_%i/result_values.%i', 4],
            'Ferr': ['det_%i/result_values.%i', 5],
            'I': ['det_%i/result_values.%i', 6],
            'I_c': ['det_%i/result_values.%i', 7],
            'I_r': ['det_%i/result_values.%i', 8],
            'Ibgr': ['det_%i/result_values.%i', 9],
            'Ibgr_c': ['det_%i/result_values.%i', 10],
            'Ibgr_r': ['det_%i/result_values.%i', 11],
            'Ierr': ['det_%i/result_values.%i', 12],
            'Ierr_c': ['det_%i/result_values.%i', 13],
            'Ierr_r': ['det_%i/result_values.%i', 14]}

DET_ATT_KEYS = ['data',
                'name']

VERSIONED_KEYS = []
                    
##############################################################################
class HdfDataFile:
    """
    Container for data stored in HDF files
    """
    def __init__(self,fname):
        #
        self.fname = fname
        self.file = h5py.File(self.fname,'r+')
        self.all_items = self.file.items()
        #self.all_items = None
        self.point = 0
        self.point_dict = {}
        self.version = 1
        
    def __del__(self):
        """
        If the object is going to be destroyed, make sure
        the current point is written back to the file first.
        
        """
        
        if self.point != 0 and self.point_dict != {}:
            try:
                self.write_point(self.point_dict, self.point)
            except:
                pass
        
        self.point = 0
        self.point_dict = {}
        try:
            self.file.flush()
            self.file.close()
        except:
            pass
        del self.file
        del self.all_items
        del self
    
    def __getitem__(self,arg):
        """
        Since usage will be hdf_object[point][arg], this reads in
        the point and returns the dictionary.
        
        """
        '''print "**arguments=", arg
        self._check_file()
        if type(arg) == types.StringType:
            #result = DEFAULT_POINT_DATA[arg]
            #str = "/pointdata/%4d/%s"  % (self.point,arg)
            str = "/%s" % arg
            #result = self.file.get(str)
            result = self.file[str]
            
            #if result == hdf_group:
            #    result = {}
            #    tunnel into group building 
            return result
        # if the file doesnt have given attribute / data field
        # how do we return error???'''
        if arg == self.point:
            return self.point_dict
        if self.point != 0 and self.point_dict != {}:
            self.write_point(self.point_dict, self.point)
        self.read_point(arg)        
        return self.point_dict
    
    def close(self):
        """
        If the object is going to be closed, make sure
        the current point is written back to the file first.
        
        """
        
        if self.point != 0 and self.point_dict != {}:
            self.write_point(self.point_dict, self.point)
        
        self.point = 0
        self.point_dict = {}
        
        self.file.flush()
        self.file.close()
    
    def delete(self, item):
        """Delete a point from the file."""
        
        del self.file[item]
    
    def get(self, num, default=None):
        """
        Acts like a dictionary's get: if num exists, returns
        the associated dictionary, otherwise returns default.
        
        """
        try:
            return self.__getitem__(num)
        except KeyError:
            return default
        except:
            raise
    
    def get_all(self, key, points=None):
        """
        Gets the value of key for every point in points.
        If points is None, gets the value for every point
        in the file. To tunnel, pass a tuple to key,
        eg HdfObject.get_all(('det_0', 'image_data'))
        
        To ensure the returned values are up to date, replaces
        the point's value (if appropriate) with the value in
        the current dictionary.
        
        """
        
        all_results = {}

        #if self.point != 0 and self.point_dict != {}:
        #    self.write_point(self.point_dict, self.point)

        if points == None:
            points = []
            for item in self.all_items:
                points.append(item[0])
        #for point in points:
        if isinstance(key, basestring):
            if key in GEN_KEYS:
                key_loc = GEN_KEYS[key]
                for point in points:
                    all_results[point] = \
                            self.file[point][key_loc[0]][key_loc[1]]
                if self.point in points:
                    all_results[self.point] = self.point_dict[key]
            elif key in ATT_KEYS:
                if key.startswith('hist'):
                    key = key % self.version
                for point in points:
                    all_results[point] = self.file[point].attrs[key]
                if self.point in points:
                    all_results[self.point] = self.point_dict[key]
            elif key in MISC_KEYS:
                for point in points:
                    all_results[point] = self.file[point][key]
                if self.point in points:
                    all_results[self.point] = self.point_dict[key]
            else:
                for point in points:
                    if key in self.file[point]['position_labels']:
                        key_loc = \
                          list(self.file[point]['position_labels']).index(key)
                        all_results[point] = \
                                self.file[point]['position_values'][key_loc]
                    elif key in self.file[point]['scaler_labels']:
                        key_loc = \
                          list(self.file[point]['scaler_labels']).index(key)
                        all_results[point] = \
                                self.file[point]['scaler_values'][key_loc]
                    else:
                        print 'Error (unrecognized key): ' , key
                if self.point in points and key in self.point_dict.keys():
                    all_results[self.point] = self.point_dict[key]
        elif isinstance(key, tuple):
            det_name = key[0]
            key = key[1]
            if key in DET_KEYS:
                key_loc = DET_KEYS[key]
                key_loc_path = key_loc[0].split('/')[1] % self.version
                for point in points:
                    all_results[point] = \
                            self.file[point][det_name][key_loc_path][key_loc[1]]
                if self.point in points:
                    all_results[self.point] = self.point_dict[det_name][key]
            elif key in DET_ATT_KEYS:
                for point in points:
                    all_results[point] = self.file[point][det_name].attrs[key]
                if self.point in points:
                    all_results[self.point] = self.point_dict[det_name][key]
            elif key.startswith('image_data'):
                for point in points:
                    try:
                        all_results[point] = self.file[point][det_name][key]
                    except:
                        pass
            elif key.startswith('corrected_image'):
                for point in points:
                    try:
                        point_image = \
                          numpy.array(self.file[point][det_str]['image_data'])
                        point_mask = self.point_dict[det_str]['bad_pixel_map']
                        all_results[point] = \
                                image_data.correct_image(point_image,
                                                         point_mask)
                    except:
                        pass
                if self.point in points:
                    all_results[self.point] = self.point_dict[det_name][key]
            else:
                print 'Error: unrecognized key'
        else:
            print 'Error: unknown key type'
        return all_results
        
    def read_point(self,num):
        """
        read data from the point to self
        
        return all data as a dictionary
        
        num should be the whole serial number string,
        eg '000328'
        
        """
        #self._check_file()
        self.point = num
        self.point_dict = {}
        for key in GEN_KEYS:
            key_loc = GEN_KEYS[key]
            self.point_dict[key] = self.file[num][key_loc[0]][key_loc[1]]
        for key in self.file[num]['position_labels']:
            key_loc = list(self.file[num]['position_labels']).index(key)
            self.point_dict[key] = self.file[num]['position_values'][key_loc]
        for key in self.file[num]['scaler_labels']:
            key_loc = list(self.file[num]['scaler_labels']).index(key)
            self.point_dict[key] = self.file[num]['scaler_values'][key_loc]
        for key in ATT_KEYS:
            if key.startswith('hist'):
                key = key % self.version
            self.point_dict[key] = self.file[num].attrs[key]
        for key in MISC_KEYS:
            self.point_dict[key] = self.file[num][key]
        current_det_num = 0
        while True:
            try:
                det_str = 'det_%i' % current_det_num
                current_det = self.file[num][det_str]
                self.point_dict[det_str] = {}
                for key in DET_KEYS:
                    key_loc = DET_KEYS[key]
                    key_loc_path = key_loc[0] % (current_det_num, self.version)
                    self.point_dict[det_str][key] = \
                            self.file[num][key_loc_path][key_loc[1]]
                for key in DET_ATT_KEYS:
                    self.point_dict[det_str][key] = \
                            self.file[num][det_str].attrs[key]
                try:
                    point_image = \
                            numpy.array(self.file[num][det_str]['image_data'])
                    self.point_dict[det_str]['image_data'] = point_image
                    point_mask = self.point_dict[det_str]['bad_pixel_map']
                    corrected_image = image_data.correct_image(point_image,
                                                               point_mask)
                    self.point_dict[det_str]['corrected_image'] = \
                                                                corrected_image
                except:
                    pass
                current_det_num += 1
            except:
                break
    
    '''def set(self, num, key, value):
        """Overwrite key in num with value."""
        pass'''
    
    def set_all(self, key, value, points=None):
        """
        Sets the value of key for every point in points.
        If points is None, sets the value for every point
        in the file. To tunnel, pass a tuple to key,
        eg HdfObject.set_all(('det_0', 'bad_pixel_map'))
        
        To ensure the set values aren't overwritten by the
        current dictionary when it's written, sets that value
        as well (if appropriate).
        
        """
        
        #if self.point != 0 and self.point_dict != {}:
        #    self.write_point(self.point_dict, self.point)
        
        if points == None:
            points = []
            for item in self.all_items:
                points.append(item[0])
        #for point in points:
        if isinstance(key, basestring):
            if key in GEN_KEYS:
                if self.point in points:
                    self.point_dict[key] = value
                key_loc = GEN_KEYS[key]
                for point in points:
                    self.file[point][key_loc[0]][key_loc[1]] = value
            elif key in self.file[num]['position_labels']:
                if self.point in points:
                    self.point_dict[key] = value
                key_loc = list(self.file[num]['position_labels']).index(key)
                for point in points:
                    self.file[point]['position_values'][key_loc] = value
            elif key in self.file[num]['scaler_labels']:
                if self.point in points:
                    self.point_dict[key] = value
                key_loc = list(self.file[num]['scaler_labels']).index(key)
                for point in points:
                    self.file[point]['scaler_values'][key_loc] = value
            elif key in ATT_KEYS:
                if self.point in points:
                    self.point_dict[key] = value
                if key.startswith('hist'):
                    key = key % self.version
                for point in points:
                    self.file[point].attrs[key] = value
            elif key in MISC_KEYS:
                if self.point in points:
                    self.point_dict[key] = value
                for point in points:
                    self.file[point][key] = value
            else:
                print 'Error: unrecognized key'
        elif isinstance(key, tuple):
            det_name = key[0]
            key = key[1]
            if key in DET_KEYS:
                if self.point in points:
                    self.point_dict[det_name][key] = value
                key_loc = DET_KEYS[key]
                key_loc_path = key_loc[0].split('/')[1] % self.version
                for point in points:
                    try:
                        self.file[point][det_name][key_loc_path][key_loc[1]] = \
                                                                           value
                    except IOError:
                        self.file[point][det_name][key_loc_path][key_loc[1]] = \
                                                              numpy.float(value)
            elif key in DET_ATT_KEYS:
                if self.point in points:
                    self.point_dict[det_name][key] = value
                for point in points:
                    self.file[point][det_name].attrs[key] = value
            elif key.startswith('image_data'):
                print 'Are you sure you want to overwrite the image data?'
                print 'If so, go into the hdf_data.py file and uncomment ' + \
                        'the lines following this message.'
                '''
                if self.point in points:
                    self.point_dict[det_name][key] = value
                for point in points:
                    try:
                        self.file[point][det_name][key] = value
                    except:
                        pass
                '''
            elif key.startswith('corrected_image'):
                pass
            else:
                print 'Error: unrecognized key'
        else:
            print 'Error: unknown key type'
        
    def version_point(self,data={}):
        """
        edit data for a given point, making a new version
        
        data is dictionary, just with the new/updated stuff
        """
        pass
    
    def write_point(self,data,num=None):
        """
        write data to file 
        
        data is a dictionary 
        """
        #self._check_file()
        #
        if num is None:
            num = self.point
        for key in data:
            if key.startswith('hist'):
                self.file[num].attrs[key] = data[key]
            elif key.startswith('det_'):
                det_dict = data[key]
                for det_key in det_dict:
                    try:
                        key_loc = DET_KEYS[det_key]
                        key_loc_path = key_loc[0].split('/')[1] % self.version
                        try:
                            self.file[num][key][key_loc_path][key_loc[1]] = \
                                                              data[key][det_key]
                        except IOError:
                            self.file[num][key][key_loc_path][key_loc[1]] = \
                                                 numpy.float(data[key][det_key])
                    except KeyError:
                        pass
            else:
                pass

##############################################################################
if __name__ == "__main__":
    #file = h5py.File('bob.h5','w')
    #file = h5py.File('C:\\Users\\biwer\\Desktop\\HDFFiles\\ProjectFULLALL.h5','w')
    #test_grp = file.create_group('test')
    #test_grp.create_dataset('ones',data=numpy.ones(10))
    #test_grp['xx'] = 'a string'
    #file.close()
    #
    #d = HdfDataFile('bob.h5')
    #testObject = HdfDataFile('C:\\Users\\biwer\\Desktop\\HDFFiles\\ProjectFULLALLv3.h5')
    testObject = HdfDataFile('C:\\Users\\biwer\\Desktop\\HDFFiles\\WithLambda.h5')
    #testObject = HdfDataFile('C:\\Users\\biwer\\Desktop\\HDFFiles\\ImageTest.h5')
    #print d['detector_1']['name']
    #print d['test']['xx']
    #x = d['test']['ones']
    #print x.value
    #for l in dir(x): print l
    #
    #d.file.close()
    
