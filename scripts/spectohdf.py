'''
Spec to HDF5 converter / parser
Author: Craig Biwer (cbiwer@uchicago.edu)
Last modified: 7/16/2012
'''

import array
import errno
import math
import os
import sys
import time

import h5py
import numpy as num
from PIL import Image

# This runs through a specfile and grabs all the data, sorted by scan.
def summarize(lines):

        # This list corresponds to the current set of parameters
        # written by spec to lines beginning #G. If that changes,
        # at least this list will need to be modified, if not more
        # code below.
        g_labs = ['g_prefer', 'g_sect', 'g_frz', 'g_haz', 'g_kaz', 'g_laz',
                    'g_zh0', 'g_zk0', 'g_z10', 'g_zh1', 'g_zk1', 'g_zl1',
                    'g_kappa', 'g_13', 'g_14', 'g_sigtau', 'g_mode1',
                    'g_mode2', 'g_mode3', 'g_mode4', 'g_mode5', 'g_21',
                    'g_aa', 'g_bb', 'g_cc', 'g_al', 'g_be', 'g_ga', 'g_aa_s',
                    'g_bb_s', 'g_cc_s', 'g_al_s', 'g_be_s', 'g_ga_s', 'g_h0',
                    'g_k0', 'g_l0', 'g_h1', 'g_k1', 'g_l1', 'g_u00', 'g_u01',
                    'g_u02', 'g_u03', 'g_u04', 'g_u05', 'g_u10', 'g_u11',
                    'g_u12', 'g_u13', 'g_u14', 'g_u15', 'g_lambda0',
                    'g_lambda1', 'g_54', 'g_55', 'g_56', 'g_57', 'g_58',
                    'g_59', 'g_60', 'g_61', 'g_62', 'g_H', 'g_K', 'g_L',
                    'g_LAMBDA', 'g_ALPHA', 'g_BETA', 'g_OMEGA', 'g_TTH',
                    'g_PSI', 'g_TAU', 'g_QAZ', 'g_NAZ', 'g_SIGMA_AZ',
                    'g_TAU_AZ', 'g_F_ALPHA', 'g_F_BETA', 'g_F_OMEGA',
                    'g_F_PSI', 'g_F_NAZ', 'g_F_QAZ', 'g_F_DEL', 'g_F_ETA',
                    'g_F_CHI', 'g_F_PHI', 'g_F_NU', 'g_F_MU', 'g_F_CHI_Z',
                    'g_F_PHI_Z', 'CUT_DEL', 'CUT_ETA', 'CUT_CHI', 'CUT_PHI',
                    'CUT_NU', 'CUT_MU', 'CUT_KETA', 'CUT_KAP', 'CUT_KPHI',
                    'g_100', 'g_101', 'g_102', 'g_103', 'g_104', 'g_105',
                    'g_106', 'g_107', 'g_108', 'g_109', 'g_110', 'g_111']

        summary = []
        lineno = 0
        (spec_name, epoch, mnames, cmnd, date, xtime, g_vals, q, p_vals, atten,
        energy, lab, aborted) = \
            (None, None, None, None, None, None, None, None, None, None,
            None, None, False)
        point_data = []
        (index, ncols, n_sline) = (0, 0, 0)
        for i in lines:
            lineno = lineno + 1
            i = i[:-1]
            # get the name of the specfile, should only appear once
            if (i.startswith('#F')):
                spec_name = i[3:]
                spec_name = os.path.split(spec_name)[-1]
            # get the starting epoch
            if (i.startswith('#E ')):
                epoch = int(i[3:])
            # get motor names: they should be at the top of the file
            # but they can be reset anywhere in the file
            elif (i.startswith('#O')):
                if i[2] == '0': mnames = ''
                mnames = mnames + i[3:]
            # get scan number
            elif (i.startswith('#S ')):
                v = i[3:].split()
                index = int(v[0])
                cmnd = i[4 + len(v[0]):]
                n_sline = lineno
            # Get the date
            elif (i.startswith('#D ')):
                date = i[3:]
            # Get some value in seconds (unclear what it is)
            elif (i.startswith('#T ')):
                xtime = i[3:]
            # Get the G values
            elif (i.startswith('#G')):
                if i[2] == '0': g_vals = ''
                g_vals = g_vals + i[3:]
            # Get Q
            # NOTE: This value is useless: it shows what Q was at the end
            # of the previous scan, not what it is going to be.
            elif (i.startswith('#Q ')):
                q = i[3:]
            # Get the motor values
            elif (i.startswith('#P')):
                if i[2] == '0': p_vals = ''
                p_vals = p_vals + i[3:]
            # Get the number of data columns
            elif (i.startswith('#N ')):
                ncols = int(i[3:])
            # Get the attenuator positions
            # NOTE: This will also match #ATTEN_FACTORS
            # If that starts getting used again, we can either change
            # this to i.startswith('#ATTEN ') or put it after an 
            # i.startswith('#ATTEN_F') check.
            elif (i.startswith('#AT')):
                atten = i[6:].strip()
            # Get the beam energy
            elif (i.startswith('#EN')):
                energy = float(i[8:])
            # Get the column header info, as well as the following data
            # Also check for comments following the data containing
            # the word 'aborted'
            elif (i.startswith('#L ')):
                lab = i[3:].split()
                ## count how many lines of 'data' we have
                ## and see if the scan was aborted
                xx = lines[lineno:]
                nl_dat = 0
                aborted = False
                for ii in xx:
                    if (ii.startswith('#S ')):
                        break
                    elif (ii.startswith('#')):
                        if ii.find('aborted') > -1:
                            aborted = True
                    elif (len(ii) > 3):
                        nl_dat = nl_dat + 1
                        point_data.append(map(float, ii.split()))
                ## append all the info...
                L_pos = -1
                L_start = '--'
                L_stop = '--'
                try:
                    L_pos = lab.index('L')
                except ValueError:
                    pass
                if L_pos > -1 and len(point_data) > 0:
                    L_start = point_data[0][L_pos]
                    L_stop = point_data[-1][L_pos]
                current_dict = {'index':index,
                                     'spec_name':spec_name,
                                     'init_epoch':epoch,
                                     'nl_start':n_sline,
                                     'cmd':cmnd,
                                     'date':date,
                                     'time':xtime,
                                     'mnames':mnames.split(),
                                     'P':map(float, p_vals.split()),
                                     'g_labs':g_labs,
                                     'G':map(float, g_vals.split()),
                                     'Q':q,
                                     'ncols':ncols,
                                     'labels':lab,
                                     'atten':atten,
                                     'energy':energy,
                                     'lineno':lineno,
                                     'aborted':aborted,
                                     'point_data':point_data,
                                     'real_L_start':L_start,
                                     'real_L_stop':L_stop,
                                     'nl_dat':nl_dat}
                summary.append(current_dict)
                (cmnd, date, xtime, g_vals, q, p_vals, atten, energy, lab,
                aborted) = \
                  (None, None, None, None, None, None, None, None, None,
                  False)
                point_data = []
                (index, ncols, n_sline) = (0, 0, 0)

        min_scan = summary[0]['index']
        max_scan = summary[0]['index']
        for i in summary:
            k = i['index']
            if (k > max_scan): max_scan = k
            if (k < min_scan): min_scan = k
            
        return summary
    
def read_image(file):
    try:
        im = Image.open(file)
        # While numpy free, this seems to only be able to handle 1D arrays:
        #arr = array.array('L', im.tostring())
        # So we'll still have to use numpy for this
        arr = num.fromstring(im.tostring(), dtype='int32')
        arr.shape = (im.size[1], im.size[0])
        return arr
    except:
        print "Error reading file: %s" % file
        return None

# Main conversion function, takes user input and converts
# specified .spc file to .mh5
def spec_to_hdf(args):
    choice = 'a'
    if len(args) == 0:
        print 'Current directory: ', os.getcwd()
        print 'Please enter the path to the specfile:'
        input = raw_input('> ')
        if not os.path.isfile(input):
            print 'Error: file not found'
            return
        print 'Please enter the name of the master file:'
        output = raw_input('> ')
        if output == '':
            output = input
        if not output.endswith('.mh5'):
            output = output + '.mh5'
        if os.path.isfile(output):
            choice = raw_input('File already exists: (A)ppend, over(W)rite, '+\
                               'or (C)ancel (A/W/C)? ').lower()
            if choice == 'c':
                return
            elif choice == '':
                choice = 'a'
    elif len(args) == 1:
        input = args[0]
        print 'Please enter the name of the master file:'
        output = raw_input('> ')
        if output == '':
            output = input
        if not output.endswith('.mh5'):
            output = output + '.mh5'
        if os.path.isfile(output):
            choice = raw_input('File already exists: (A)ppend, over(W)rite, '+\
                               'or (C)ancel (A/W/C)? ').lower()
            if choice == 'c':
                return
            elif choice == '':
                choice = 'a'
    elif len(args) == 2:
        input = args[0]
        if not os.path.isfile(input):
            print 'Error: file not found'
            return
        output = args[1]
        if output == '':
            output = input
        if not output.endswith('.mh5'):
            output = output + '.mh5'
        if os.path.isfile(output):
            choice = raw_input('File already exists: (A)ppend, over(W)rite, '+\
                               'or (C)ancel (A/W/C)? ').lower()
            if choice == 'c':
                return
            elif choice == '':
                choice = 'a'
    else:
        options = args[2:]
        choice = None
        if '-a' in options:
            choice = 'a'
        elif '-w' in options:
            choice = 'w'
        input = args[0]
        if not os.path.isfile(input):
            print 'Error: file not found'
            return
        output = args[1]
        if output == '':
            output = input
        if not output.endswith('.mh5'):
            output = output + '.mh5'
        if os.path.isfile(output):
            if choice == 'c':
                return
    
    time1 = time.time()
    
    spec_dir = os.path.split(input)[0]
    if spec_dir == '':
        spec_dir = '.'
    spec_name = os.path.split(input)[-1]
    image_dir = spec_dir + '\\images\\%s\\' % spec_name[:-4]
    this_file = open(input)
    lines = this_file.readlines()
    this_file.close()
    summary = summarize(lines)

    # Obtain a lock for the HDF file. If the file is already
    # locked, retry for 10 seconds before returning.
    lockFile = FileLock(output, timeout=1)
    try:
        print 'Attempting to lock file...'
        lockFile.acquire()
        print 'Lock acquired'
    except FileLockException as e:
        print 'Error: ' + str(e)
        return
    
    master_file = h5py.File(output, choice)
    spec_group = master_file.require_group(spec_name)
    
    for scan in summary:
        scan_group = spec_group.require_group(str(scan['index']))
        try:
            if len(scan_group['point_data']) == scan['nl_dat']:
                print 'Skipping scan' + str(scan['index'])
                continue
            del scan_group['point_data']
        except KeyError:
            pass
        scan_group.create_dataset('point_data', data=scan['point_data'])
        try:
            del scan_group['point_labs']
        except KeyError:
            pass
        scan_group.create_dataset('point_labs', data=scan['labels'])
        try:
            del scan_group['param_labs']
        except KeyError:
            pass
        scan_group.create_dataset('param_labs',
                                 data=scan['g_labs'] + scan['mnames'])
        try:
            del scan_group['param_data']
        except KeyError:
            pass
        scan_group.create_dataset('param_data', data=scan['G'] + scan['P'])
        # Scan types: ['a2scan', 'a4scan', 'ascan', 'Escan', 'hkcircle',
        #              'hklmesh', 'hklrock', 'hklscan', 'loopscan', 'mesh',
        #              'powder_scan', 'rodscan', 'rodscan_interp',
        #              'timescan', 'xascan']
        # Currently included: a2scan, a4scan, ascan, Escan, hklscan, rodscan,
        #                     timescan
        # Parses according to the command given
        if scan.get('cmd', '').split()[0] == 'a2scan':
            for key in scan.keys():
                if key == 'cmd':
                    scan_group.attrs['cmd'] = scan[key]
                    split_cmd = scan[key].split()
                    scan_group.attrs['s_type'] = split_cmd[0]
                    scan_group.attrs['motor1'] = split_cmd[1]
                    scan_group.attrs['m1_start'] = float(split_cmd[2])
                    scan_group.attrs['m1_stop'] = float(split_cmd[3])
                    scan_group.attrs['motor2'] = split_cmd[4]
                    scan_group.attrs['m2_start'] = float(split_cmd[5])
                    scan_group.attrs['m2_stop'] = float(split_cmd[6])
                    scan_group.attrs['num_points'] = float(split_cmd[7])
                    scan_group.attrs['count_time'] = float(split_cmd[8])
                elif key not in ['labels', 'point_data',
                               'g_labs', 'mnames', 'G', 'P']:
                    scan_key = scan[key]
                    if scan_key is None:
                        scan_key = 'None'
                    scan_group.attrs[key] = scan_key
                    '''if key.startswith('date'):
                        try:
                            scan_group.attrs['init_epoch'] = \
                                        time.mktime(time.strptime(scan_key))
                        except:
                            pass'''
        elif scan.get('cmd', '').split()[0] == 'a4scan':
            for key in scan.keys():
                if key == 'cmd':
                    scan_group.attrs['cmd'] = scan[key]
                    split_cmd = scan[key].split()
                    scan_group.attrs['s_type'] = split_cmd[0]
                    scan_group.attrs['motor1'] = split_cmd[1]
                    scan_group.attrs['m1_start'] = float(split_cmd[2])
                    scan_group.attrs['m1_stop'] = float(split_cmd[3])
                    scan_group.attrs['motor2'] = split_cmd[4]
                    scan_group.attrs['m2_start'] = float(split_cmd[5])
                    scan_group.attrs['m2_stop'] = float(split_cmd[6])
                    scan_group.attrs['motor3'] = split_cmd[7]
                    scan_group.attrs['m3_start'] = float(split_cmd[8])
                    scan_group.attrs['m3_stop'] = float(split_cmd[9])
                    scan_group.attrs['motor4'] = split_cmd[10]
                    scan_group.attrs['m4_start'] = float(split_cmd[11])
                    scan_group.attrs['m4_stop'] = float(split_cmd[12])
                    scan_group.attrs['num_points'] = float(split_cmd[13])
                    scan_group.attrs['count_time'] = float(split_cmd[14])
                elif key not in ['labels', 'point_data',
                               'g_labs', 'mnames', 'G', 'P']:
                    scan_key = scan[key]
                    if scan_key is None:
                        scan_key = 'None'
                    scan_group.attrs[key] = scan_key
                    '''if key.startswith('date'):
                        try:
                            scan_group.attrs['epoch'] = \
                                        time.mktime(time.strptime(scan_key))
                        except:
                            pass'''
        elif scan.get('cmd', '').split()[0] == 'ascan':
            for key in scan.keys():
                if key == 'cmd':
                    scan_group.attrs['cmd'] = scan[key]
                    split_cmd = scan[key].split()
                    scan_group.attrs['s_type'] = split_cmd[0]
                    scan_group.attrs['motor1'] = split_cmd[1]
                    scan_group.attrs['m1_start'] = float(split_cmd[2])
                    scan_group.attrs['m1_stop'] = float(split_cmd[3])
                    scan_group.attrs['num_points'] = float(split_cmd[4])
                    scan_group.attrs['count_time'] = float(split_cmd[5])
                elif key not in ['labels', 'point_data',
                               'g_labs', 'mnames', 'G', 'P']:
                    scan_key = scan[key]
                    if scan_key is None:
                        scan_key = 'None'
                    scan_group.attrs[key] = scan_key
                    '''if key.startswith('date'):
                        try:
                            scan_group.attrs['epoch'] = \
                                        time.mktime(time.strptime(scan_key))
                        except:
                            pass'''
        elif scan.get('cmd', '').split()[0] == 'Escan':
            for key in scan.keys():
                if key == 'cmd':
                    scan_group.attrs['cmd'] = scan[key]
                    split_cmd = scan[key].split()
                    scan_group.attrs['s_type'] = split_cmd[0]
                    scan_group.attrs['energy_start'] = float(split_cmd[1])
                    scan_group.attrs['energy_stop'] = float(split_cmd[2])
                    scan_group.attrs['num_points'] = float(split_cmd[3])
                    scan_group.attrs['count_time'] = float(split_cmd[4])
                elif key.startswith('Q'):
                    scan_group.attrs['h_val'] = \
                                        round(float(scan[key].split()[0]), 1)
                    scan_group.attrs['k_val'] = \
                                        round(float(scan[key].split()[1]), 1)
                elif key not in ['labels', 'point_data',
                               'g_labs', 'mnames', 'G', 'P', 'Q']:
                    scan_key = scan[key]
                    if scan_key is None:
                        scan_key = 'None'
                    scan_group.attrs[key] = scan_key
                    '''if key.startswith('date'):
                        try:
                            scan_group.attrs['epoch'] = \
                                        time.mktime(time.strptime(scan_key))
                        except:
                            pass'''
        elif scan.get('cmd', '').split()[0] == 'hklscan':
            for key in scan.keys():
                if key == 'cmd':
                    scan_group.attrs['cmd'] = scan[key]
                    split_cmd = scan[key].split()
                    scan_group.attrs['s_type'] = split_cmd[0]
                    scan_group.attrs['h_val'] = round(float(split_cmd[1]), 1)
                    scan_group.attrs['h_start'] = float(split_cmd[1])
                    scan_group.attrs['h_stop'] = float(split_cmd[2])
                    scan_group.attrs['k_val'] = round(float(split_cmd[3]), 1)
                    scan_group.attrs['k_start'] = float(split_cmd[3])
                    scan_group.attrs['k_stop'] = float(split_cmd[4])
                    scan_group.attrs['L_start'] = float(split_cmd[5])
                    scan_group.attrs['L_stop'] = float(split_cmd[6])
                    scan_group.attrs['num_points'] = float(split_cmd[7])
                    scan_group.attrs['count_time'] = float(split_cmd[8])
                elif key not in ['labels', 'point_data',
                               'g_labs', 'mnames', 'G', 'P']:
                    scan_key = scan[key]
                    if scan_key is None:
                        scan_key = 'None'
                    scan_group.attrs[key] = scan_key
                    '''if key.startswith('date'):
                        try:
                            scan_group.attrs['epoch'] = \
                                        time.mktime(time.strptime(scan_key))
                        except:
                            pass'''
        elif scan.get('cmd', '').split()[0] == 'rodscan':
            for key in scan.keys():
                if key == 'cmd':
                    scan_group.attrs['cmd'] = scan[key]
                    split_cmd = scan[key].split()
                    scan_group.attrs['s_type'] = split_cmd[0]
                    scan_group.attrs['h_val'] = float(split_cmd[1])
                    scan_group.attrs['k_val'] = float(split_cmd[2])
                    scan_group.attrs['L_start'] = float(split_cmd[3])
                    scan_group.attrs['L_stop'] = float(split_cmd[4])
                    scan_group.attrs['L_space'] = float(split_cmd[5])
                    scan_group.attrs['max_time'] = float(split_cmd[6])
                    scan_group.attrs['peak_pos'] = float(split_cmd[7])
                    scan_group.attrs['peak_space'] = float(split_cmd[8])
                elif key not in ['labels', 'point_data',
                               'g_labs', 'mnames', 'G', 'P']:
                    scan_key = scan[key]
                    if scan_key is None:
                        scan_key = 'None'
                    scan_group.attrs[key] = scan_key
                    '''if key.startswith('date'):
                        try:
                            scan_group.attrs['epoch'] = \
                                        time.mktime(time.strptime(scan_key))
                        except:
                            pass'''
            # Calculate the HK distance
            h_val = scan_group.attrs['h_val']
            k_val = scan_group.attrs['k_val']
            h_val_d = list(scan_group['param_labs']).index('g_aa_s')
            k_val_d = list(scan_group['param_labs']).index('g_bb_s')
            t_val_d = list(scan_group['param_labs']).index('g_ga_s')
            h_val_d = float(list(scan_group['param_data'])[h_val_d])
            k_val_d = float(list(scan_group['param_data'])[k_val_d])
            t_val_d = float(list(scan_group['param_data'])[t_val_d])
            t_val_d = math.radians(180 - t_val_d)
            hk_dist = round(((h_val * h_val_d)**2 + \
                             (k_val * k_val_d)**2 - \
                              2 * h_val * h_val_d * k_val * k_val_d * \
                              math.cos(t_val_d))**0.5, 8)
            scan_group.attrs['hk_dist'] = hk_dist
        elif scan.get('cmd', '').split()[0] == 'timescan':
            for key in scan.keys():
                if key == 'cmd':
                    scan_group.attrs['cmd'] = scan[key]
                    split_cmd = scan[key].split()
                    scan_group.attrs['s_type'] = split_cmd[0]
                    scan_group.attrs['count_time'] = float(split_cmd[1])
                    scan_group.attrs['time_space'] = float(split_cmd[2])
                elif key not in ['labels', 'point_data',
                               'g_labs', 'mnames', 'G', 'P']:
                    scan_key = scan[key]
                    if scan_key is None:
                        scan_key = 'None'
                    scan_group.attrs[key] = scan_key
                    '''if key.startswith('date'):
                        try:
                            scan_group.attrs['epoch'] = \
                                        time.mktime(time.strptime(scan_key))
                        except:
                            pass'''
        else:
            for key in scan.keys():
                if key not in ['labels', 'point_data',
                               'g_labs', 'mnames', 'G', 'P']:
                    scan_key = scan[key]
                    if scan_key is None:
                        scan_key = 'None'
                    scan_group.attrs[key] = scan_key
                    '''if key.startswith('date'):
                        try:
                            scan_group.attrs['epoch'] = \
                                        time.mktime(time.strptime(scan_key))
                        except:
                            pass'''
        # Set the image directory for the scan
        this_dir = image_dir + 'S%03d\\' % scan['index']
        # Print the directory (to give users a sense of progress made)
        print this_dir
        if os.path.isdir(this_dir):
            image_data = []
            for image_file in os.listdir(this_dir):
                try:
                    image_value = None
                    image_path = os.path.join(this_dir, image_file)
                    if image_path.endswith('.tif'):
                        image_value = read_image(image_path)
                        if image_value != None:
                            image_data.append(image_value)
                        else:
                            holding = image_data[-1].copy()
                            holding.fill(-1)
                            image_data.append(holding)
                except:
                    print 'Error reading image ' + image_path
                    raise
            if image_data != []:
                try:
                    del scan_group['image_data']
                except KeyError:
                    pass
                scan_group.create_dataset('image_data', data=image_data,
                                         compression='szip')
    
    master_file.close()
    
    lockFile.release()
    print 'Lock released'
    
    time2 = time.time()
    print 'Total time:', (time2-time1)/60, 'minutes.'

"""
File Locker
Author: Evan Fosmark
http://www.evanfosmark.com/2009/01/cross-platform-file-locking-support-in-python
Last modified: 7.16.2012 by Craig Biwer (cbiwer@uchicago.edu)
"""
 
class FileLockException(Exception):
    pass
 
class FileLock(object):
    """ A file locking mechanism that has context-manager support so 
        you can use it in a with statement. This should be relatively cross
        compatible as it doesn't rely on msvcrt or fcntl for the locking.
    """
 
    def __init__(self, file_name, timeout=10, delay=.05):
        """ Prepare the file locker. Specify the file to lock and optionally
            the maximum timeout and the delay between each attempt to lock.
        """
        self.is_locked = False
        #self.lockfile = os.path.join(os.getcwd(), "%s.lock" % file_name)
        self.lockfile = file_name + '.lock'
        self.file_name = file_name
        self.timeout = timeout
        self.delay = delay
 
 
    def acquire(self):
        """ Acquire the lock, if possible. If the lock is in use, it check again
            every `wait` seconds. It does this until it either gets the lock or
            exceeds `timeout` number of seconds, in which case it throws 
            an exception.
        """
        start_time = time.time()
        while True:
            try:
                self.fd = os.open(self.lockfile, os.O_CREAT|os.O_EXCL|os.O_RDWR)
                os.write(self.fd, os.environ['USERNAME'] + '\n')
                os.write(self.fd, os.environ['COMPUTERNAME'] + '\n')
                os.write(self.fd, time.ctime(time.time()))
                break;
            except OSError as e:
                if e.errno != errno.EEXIST and e.errno != errno.EACCES:
                    raise 
                if (time.time() - start_time) >= self.timeout:
                    if e.errno == errno.EEXIST:
                        raise FileLockException("Timeout occured.")
                    else:
                        raise FileLockException("Access denied.")
                time.sleep(self.delay)
        self.is_locked = True
 
 
    def release(self):
        """ Get rid of the lock by deleting the lockfile. 
            When working in a `with` statement, this gets automatically 
            called at the end.
        """
        if self.is_locked:
            os.close(self.fd)
            os.unlink(self.lockfile)
            self.is_locked = False
 
 
    def __enter__(self):
        """ Activated when used in the with statement. 
            Should automatically acquire a lock to be used in the with block.
        """
        if not self.is_locked:
            self.acquire()
        return self
 
 
    def __exit__(self, type, value, traceback):
        """ Activated at the end of the with statement.
            It automatically releases the lock if it isn't locked.
        """
        if self.is_locked:
            self.release()
 
 
    def __del__(self):
        """ Make sure that the FileLock instance doesn't leave a lockfile
            lying around.
        """
        self.release()

if __name__ == '__main__':    
    args = sys.argv[1:]
    spec_to_hdf(args)
    
