#
# Read/Write CARS mca/med files
#

import Mca
import Med
import Numeric
import string
import os

#########################################################################
# Read Files
#########################################################################

##########################################################################
def read_med(file, netcdf=0):
    """
    Reads a disk file into an Med object. The file contains the information
    from the Med object which it makes sense to store permanently, but does
    not contain all of the internal state information for the Med.

    Inputs:
        file: The name of the disk file to read.
    """
    if (netcdf != 0):
        #r = Mca.read_netcdf_file(file)
        r = read_netcdf_file(file)
    else:
        #r = Mca.read_ascii_file(file)
        r = read_ascii_file(file)

    if r == None: return None
    
    n_detectors = r['n_detectors']
    path, fname = os.path.split(file)
    med = Med.Med(n_detectors=n_detectors, name=fname)

    med.mcas = []
    for i in range(n_detectors):
        med.mcas.append(Mca.Mca())
        med.mcas[i].set_rois(r['rois'][i])
        med.mcas[i].set_data(r['data'][i])
        med.mcas[i].set_name(fname + ':' + str(i))
    med.set_elapsed(r['elapsed'])
    med.set_calibration(r['calibration'])
    med.set_environment(r['environment'])

    return med

#############################################################################
def read_mca(file, netcdf=0, detector=0):
    """
    Reads a disk file into an MCA object.  If the netcdf=1 flag is set it
    reads a netcdf file, else it assumes the file is ASCII.
    If the data file has multiple detectors then the detector keyword can be
    used to specify which detector data to return.

    Inputs:
     file: The name of the disk file to read.
        
    Keywords:
     netcdf: Set this flag to read files written in netCDF format, otherwise
             the routine assumes that the file is in ASCII format.
             See the documentation for Mca.write_ascii_file and
             Mca.write_netcdf_file for information on the formats.

     detector: Specifies which detector to read if the file has multiple detectors.
        
    Example:
     mca = read_mca('mca.001')
    """
    if (netcdf != 0):
        r = read_netcdf_file(file)
    else:
        r = read_ascii_file(file)

    if r == None: return None

    path, fname = os.path.split(file)
    mca = Mca.Mca(name=fname)
    mca.calibration = r['calibration'][detector]
    mca.data = r['data'][detector]
    mca.elapsed = r['elapsed'][detector]
    mca.rois = r['rois'][detector]
    mca.environment = r['environment']

    return mca

########################################################################
def read_ascii_file(file):
    """
    Reads a disk file.  The file format is a tagged ASCII format.
    The file contains the information from the Mca object which it makes sense
    to store permanently, but does not contain all of the internal state
    information for the Mca.  This procedure reads files written with
    write_ascii_file().

    Inputs:
        file: The name of the disk file to read.
            
    Outputs:
        Returns a dictionary of the following type:
        'n_detectors': int,
        'calibration': [McaCalibration()],
        'elapsed':     [McaElapsed()],
        'rois':        [[McaROI()]]
        'data':        [Numeric.array]
        'environment': [[McaEnvironment()]]
        
    Example:
        m = read_ascii_file('mca.001')
        m['elapsed'][0].real_time
    """
    try:
        fp = open(file, 'r')
    except:
        print "File not found"
        return None
    
    line = ''
    start_time = ''
    data = None
 
    environment = []
    n_detectors = 1  # Assume single element data
    elapsed = [Mca.McaElapsed()]
    calibration = [Mca.McaCalibration()]
    rois = [[]]
    while(1):
        line = fp.readline()
        if (line == ''): break
        pos = string.find(line, ' ')
        if (pos == -1): pos = len(line)
        tag = line[0:pos]
        value = string.strip(line[pos:])
        values = string.split(value)
        if (tag == 'VERSION:'):
             pass
        elif (tag == 'DATE:'):  
            start_time = value
        elif (tag == 'ELEMENTS:'):
            n_detectors  = int(value)
            for det in range(1, n_detectors):
                 elapsed.append(Mca.McaElapsed())
                 calibration.append(Mca.McaCalibration())
                 rois.append([])
        elif (tag == 'CHANNELS:'):
            nchans = int(value)
        elif (tag == 'ROIS:'):
            nrois = []
            for d in range(n_detectors):
                nrois.append(int(values[d]))
            max_rois = max(nrois)
            for d in range(n_detectors):
                for r in range(nrois[d]):
                    rois[d].append(Mca.McaROI())
        elif (tag == 'REAL_TIME:'):
            for d in range(n_detectors):
                elapsed[d].start_time = start_time
                elapsed[d].real_time = float(values[d])
        elif (tag == 'LIVE_TIME:'):  
            for d in range(n_detectors):
                elapsed[d].live_time = float(values[d])
        elif (tag == 'INPUT_COUNTS:'):  
            for d in range(n_detectors):
                elapsed[d].input_counts = float(values[d])
        elif (tag == 'CAL_OFFSET:'):
            for d in range(n_detectors):
                calibration[d].offset = float(values[d])
        elif (tag == 'CAL_SLOPE:'):
            for d in range(n_detectors):
                calibration[d].slope = float(values[d])
        elif (tag == 'CAL_QUAD:'):  
            for d in range(n_detectors):
                calibration[d].quad = float(values[d])
        elif (tag == 'TWO_THETA:'):
            for d in range(n_detectors):
                calibration[d].two_theta = float(values[d])
        elif (tag == 'ENVIRONMENT:'):
            env = Mca.McaEnvironment()
            p1 = string.find(value, '=')
            env.name = value[0:p1]
            p2 = string.find(value[p1+2:], '"')
            env.value = value[p1+2: p1+2+p2]
            env.description = value[p1+2+p2+3:-1]
            environment.append(env)
        elif (tag == 'DATA:'):
            data = []
            for d in range(n_detectors):
                data.append(Numeric.zeros(nchans, 'i'))
            for chan in range(nchans):
                line = fp.readline()
                counts = string.split(line)
                for d in range(n_detectors):
                    data[d][chan]=int(counts[d])

            for d in range(n_detectors):
                elapsed[d].total_counts = sum(data[d])

        else:
            for i in range(max_rois):
                roi = 'ROI_'+str(i)+'_'
                if (tag == roi+'LEFT:'):
                    for d in range(n_detectors):
                        if (i < nrois[d]):
                            rois[d][i].left = int(values[d])
                    break
                elif (tag == roi+'RIGHT:'):
                    for d in range(n_detectors):
                        if (i < nrois[d]):
                            rois[d][i].right = int(values[d])
                    break
                elif (tag == roi+'LABEL:'):
                    labels = string.split(value, '&')
                    for d in range(n_detectors):
                        if (i < nrois[d]):
                            rois[d][i].label = string.strip(labels[d])
                    break
            else:
                print 'Unknown tag = '+tag+' in file: ' + file + '.'

    # Make sure DATA array is defined, else this was not a valid data file
    if (data == None): print 'Not a valid data file: ' + file + '.'
    fp.close()
    # Built dictionary to return
    r = {}
    r['n_detectors'] = n_detectors
    r['calibration'] = calibration
    r['elapsed'] = elapsed
    r['rois'] = rois
    r['data'] = data
    r['environment'] = environment
    return r

#########################################################################
# Write Files
#########################################################################

#########################################################################
def write_file(detector, file, netcdf=0):
    """
    Writes Mca or Med objects to a disk file.
    
    It calls Mca.write_netcdf_file if the netcdf keyword flg is set,
    Note that users who want to read such files with Python are strongly
    encouraged to use Mca.read_file()

    Inputs:
        file: The name of the disk file to write.
            
    Keywords:
        netcdf: Set this flag to write the file in netCDF format, otherwise
                the file is written in ASCII format.  See the documentation
                for Mca.write_ascii_file and Mca.write_netcdf_file for 
                information on the formats.

    Example:
        write_file(mca,'mca.001')
    """
    # Call the get_xxx() methods to make sure things are up to date
    data        = detector.get_data()
    calibration = detector.get_calibration()
    elapsed     = detector.get_elapsed()
    # note rois should always be in channel units
    rois        = detector.get_rois()
    environment = detector.get_environment()

    if (netcdf != 0):
        write_netcdf_file(file, data, calibration, elapsed, rois, environment)
    else:
        write_ascii_file(file, data, calibration, elapsed, rois, environment)

#######################################################################
def write_ascii_file(file, data, calibration, elapsed, rois, environment):
    """
    Writes Mca or Med data to a disk file.  The file 
    format is a tagged ASCII format.  The file contains the information 
    from the Mca object which it makes sense to store permanently, but 
    does not contain all of the internal state information for the Mca.  
    Files written with this routine can be read with read_ascii_file(), which
    is called by Mca.read_file() if the netcdf flag is 0.

    This procedure is typically not called directly, but is called
    by Mca.write_file if the netcdf=1 keyword is not used.

    This function can be used for writing for Mca objects, in which case
    each input parameter is an object, such as McaElapsed, etc.
    It can also be used for writing Med objects, in which case each input
    parameter is a list.
    
    If the rank of data is 2 then this is an Med, and the number of detectors
    is the first dimension of data

    Inputs:
        file: The name of the disk file to write.
            
        data: The data to write.  Either 1-D array or list of 1-D arrays.

        calibration: An object of type McaCalibration, or a list of such objects.

        elapsed: An object of type McaElapsed, or a list of such objects.

        rois: A list of McaROI objects, or a list of lists of such objects.
        
        environment: A list of McaEnvironment objects, or a list of lists of such objects.
    """
    if (Numeric.rank(data) == 2):
        n_det = len(data)
    else:
        n_det = 1
    fformat = '%f ' * n_det
    eformat = '%e ' * n_det
    iformat = '%d ' * n_det
    sformat = '%s ' * n_det
    if (n_det == 1):
        # For convenience we convert all attributes to lists
        data = [data]
        rois = [rois]
        calibration = [calibration]
        elapsed = [elapsed]
    nchans = len(data[0])
    start_time = elapsed[0].start_time

    fp = open(file, 'w')
    fp.write('VERSION:    '+'3.1'+'\n')
    fp.write('ELEMENTS:   '+str(n_det)+'\n')
    fp.write('DATE:       '+str(start_time)+'\n')
    fp.write('CHANNELS:   '+str(nchans)+'\n')

    nrois = []
    for roi in rois:
        nrois.append(len(roi))
    fp.write('ROIS:       '+(iformat % tuple(nrois))+'\n') 
    real_time=[]; live_time=[] ; input_counts=[]
    for e in elapsed:
        real_time.append(e.real_time)
        live_time.append(e.live_time)
        input_counts.append(e.input_counts)
    fp.write('REAL_TIME:  '+(fformat % tuple(real_time))+'\n')
    fp.write('LIVE_TIME:  '+(fformat % tuple(live_time))+'\n')
    fp.write('INPUT_COUNTS: '+(fformat % tuple(input_counts))+'\n')
    offset=[]; slope=[]; quad=[]; two_theta=[]
    for c in calibration:
        offset.append(c.offset)
        slope.append(c.slope)
        quad.append(c.quad)
        two_theta.append(c.two_theta)
    fp.write('CAL_OFFSET: '+(eformat % tuple(offset))+'\n')
    fp.write('CAL_SLOPE: '+(eformat % tuple(slope))+'\n')
    fp.write('CAL_QUAD: '+(eformat % tuple(quad))+'\n')
    fp.write('TWO_THETA: '+(fformat % tuple(two_theta))+'\n')

    # note ROIS should always be in channel units!
    for i in range(max(nrois)):
        num = str(i)
        left=[]; right=[]; label=[]
        for d in range(n_det):
            if (i < nrois[d]):
                left.append(rois[d][i].left)
                right.append(rois[d][i].right)
                label.append(rois[d][i].label + '&')
            else:
                left.append(0)
                right.append(0)
                label.append(' &')
        fp.write('ROI_'+num+'_LEFT:  '+(iformat % tuple(left))+'\n')
        fp.write('ROI_'+num+'_RIGHT:  '+(iformat % tuple(right))+'\n')
        fp.write('ROI_'+num+'_LABEL:  '+(sformat % tuple(label))+'\n')
    for e in environment:
        fp.write('ENVIRONMENT: '       + str(e.name) +
                                        '="'  + str(e.value) +
                                        '" (' + str(e.description) + ')\n')
    fp.write('DATA: \n')
    counts = Numeric.zeros(n_det)
    for i in range(nchans):
        for d in range(n_det):
            counts[d]=data[d][i]
        fp.write((iformat % tuple(counts))+'\n')
    fp.close()

#########################################################################
def increment_filename(old_file):
   """
   Increments the file extension if it is numeric.  It preserves the number of
   characters in the extension.
   
   Examples:
      print increment_filename('test.001')
         test.002
      print increment_filename('test')
         test
      print increment_filename('file.1')
         file.2
   """
   dot = old_file.rfind('.')
   if (dot == -1): return old_file
   if (dot+1 == len(old_file)): return old_file

   ext = old_file[dot+1:]
   file = old_file[0:dot+1]
   nc = str(len(ext))
   try:
      ext = int(ext)+1      # Convert to number, add one, catch error
      format = '%' + nc + '.' + nc + 'd'
      ext = (format % ext)
      new_file = file + ext
      return new_file
   except:
      return old_file
