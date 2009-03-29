########################################################################
"""
Read/Write CARS mca/med files
Original written by Mark Rivers, GSECARS
See http://cars9.uchicago.edu/software/python/index.html

--------------
 Modifications
--------------
- Modified for Tdl, tpt

"""

#########################################################################

import Mca
import Med
import ROI
import numpy as Num
import string
import os

########################################################################
class Environment:
    """
    The "environment" or related parameters for a detector.  These might include
    things like motor positions, temperature, anything that describes the
    experiment.

    Fields:
        .name         # A string name of this parameter, e.g. "13IDD:m1"
        .description  # A string description of this parameter, e.g. "X stage"
        .value        # A string value of this parameter,  e.g. "14.223"
    """
    def __init__(self, name='', value='', description=''):
        self.name        = name
        self.value       = value
        self.description = description

##############################################################################
def read_med_files(file_prefix,start=0,end=100,nfmt=3,
                   bad_mca_idx=[],total=True,align=True,correct=True,tau=None):
    """
    Read multiple files given prefix, stanrt and end numbers and fmt len
    """
    med = []
    format = '%' + str(nfmt) + '.' + str(nfmt) + 'd'
    for j in range(start,end+1):
        ext   = format % j
        file  = file_prefix + '.' + ext
        tmp   = read_med(file,bad_mca_idx=bad_mca_idx,total=total,
                         align=align, correct=correct,tau=tau)
        med.append(tmp)
    return med

##############################################################################
def read_med(file,bad_mca_idx=[],total=True,align=True,correct=True,tau=None):
    """
    Reads a disk file into an Med object. The file contains the information
    from the Med object which makes sense to store permanently, but does
    not contain all of the internal state information for the Med.

    Notes
    - If the file includes "Environment data" these wil lbe included
      in the returned object as med.environment
    - If the file includes ROI's the will be included with the mca's
      in the returned object as med.mca[i].rois

    Inputs:
        file: The name of the disk file to read.
    """
    r = read_ascii_file(file)
    if r == None: return None
    
    n_detectors = r['n_detectors']
    path, fname = os.path.split(file)
    med = Med.Med(name=fname, mcas=r['mcas'], bad_mca_idx=bad_mca_idx,
                  total=total, align=align, correct=correct)

    if tau != None:
        med.update_correction(tau=tau)
        
    # below are set only if defined in file    
    if r['max_rois'] > 0:
        rois = r['rois']
        for d in range(n_detectors):
            med.mcas[d].rois = rois[d]
    if len(r['environment']) > 0:
        med.environment = r['environment']

    return med

##############################################################################
def read_mca_files(file_prefix,start=0,end=100,nfmt=3,
                   detector=0,tau=None):
    """
    Read multiple files given prefix, stanrt and end numbers and fmt len
    """
    mca = []
    format = '%' + str(nfmt) + '.' + str(nfmt) + 'd'
    for j in range(start,end+1):
        ext   = format % j
        file  = file_prefix + '.' + ext
        tmp   = read_mca(file,detector=detector,tau=tau)
        mca.append(tmp)
    return med

#############################################################################
def read_mca(file, detector=0, tau=None):
    """
    Reads a disk file into an MCA object.  
    If the data file has multiple detectors then the detector
    keyword can be used to specify which detector data to return.

    Notes.
    - If the file includes "Environment data" these wil lbe included
      in the returned object as mca.environment
    - If the file includes ROI's the will be included with the mca's
      in the returned object as mca.rois

    Inputs:
     file: The name of the disk file to read.
     detector: Index of detector to read
    """
    r = read_ascii_file(file)
    if r == None: return None

    path, fname = os.path.split(file)
    mca = r['mcas'][detector]
    mca.update_correction(tau=tau)

    # below are set only if defined in file    
    if r['max_rois'] > 0:
        rois = r['rois']
        mca.rois = rois[detector]
    if len(r['environment']) > 0:
        mca.environment = r['environment']

    return mca

########################################################################
def read_ascii_file(file):
    """
    Reads a disk file.  The file is a tagged ASCII format.
    The file contains the information from the Mca object which it makes sense
    to store permanently, but does not contain all of the internal state
    information for the Mca.  This procedure reads files written with
    write_ascii_file().

    Inputs:
        file: The name of the disk file to read.
            
    Outputs:
        Returns a dictionary of the following type:
        r['n_detectors'] = n_detectors
        r['mcas'] = [Mca.Mca()]
        r['rois'] = [Roi.Roi()]
        r['environment'] = [Environment]
        
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
    n_detectors = 0
    nchans = 0
    mcas = []
    max_rois = 0
    rois = []
    environment = []
    data = None

    while(1):
        line = fp.readline()
        if (line == ''): break
        pos = string.find(line, ' ')
        if (pos == -1): pos = len(line)
        tag = line[0:pos]
        value = string.strip(line[pos:])
        values = string.split(value)

        # scan through tags
        # NOte the elements and channels tags
        # should be at the top of the file,
        # elements comes before channels!
        if (tag == 'VERSION:'):
             pass
        elif (tag == 'DATE:'):  
            start_time = value
        elif (tag == 'ELEMENTS:'):
            n_detectors  = int(value)
        elif (tag == 'CHANNELS:'):
            nchans = int(value)
            for det in range(n_detectors):
                name = 'mca%s' % str(det)
                mcas.append(Mca.Mca(name=name,nchans=nchans))
        elif (tag == 'REAL_TIME:'):
            for d in range(n_detectors):
                mcas[d].start_time = start_time
                mcas[d].real_time = float(values[d])
        elif (tag == 'LIVE_TIME:'):  
            for d in range(n_detectors):
                mcas[d].live_time = float(values[d])
        elif (tag == 'INPUT_COUNTS:'):  
            for d in range(n_detectors):
                mcas[d].input_counts = float(values[d])
        elif (tag == 'TAU:'):  
            for d in range(n_detectors):
                mcas[d].tau = float(values[d])
        elif (tag == 'CAL_OFFSET:'):
            for d in range(n_detectors):
                mcas[d].offset = float(values[d])
        elif (tag == 'CAL_SLOPE:'):
            for d in range(n_detectors):
                mcas[d].slope = float(values[d])
        elif (tag == 'CAL_QUAD:'):  
            for d in range(n_detectors):
                mcas[d].quad = float(values[d])
        elif (tag == 'TWO_THETA:'):
            for d in range(n_detectors):
                mcas[d].two_theta = float(values[d])
        # Note 'ROIS' tag must come before the 'ROI_' tags!
        elif (tag == 'ROIS:'):
            nrois = []
            for d in range(n_detectors):
                rois.append([])
                nrois.append(int(values[d]))
            max_rois = max(nrois)
            if max_rois > 0:
                for d in range(n_detectors):
                    for r in range(nrois[d]):
                        rois[d].append(ROI.ROI())
        # parse rois
        elif (tag[0:4] == 'ROI_'):
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
        elif (tag == 'ENVIRONMENT:'):
            env = Environment()
            p1 = string.find(value, '=')
            env.name = value[0:p1]
            p2 = string.find(value[p1+2:], '"')
            env.value = value[p1+2: p1+2+p2]
            env.description = value[p1+2+p2+3:-1]
            environment.append(env)
        # DATA should be the final tag
        elif (tag == 'DATA:'):
            for chan in range(nchans):
                line = fp.readline()
                counts = string.split(line)
                for d in range(n_detectors):
                    mcas[d].data[chan]=int(counts[d])
            for d in range(n_detectors):
                mcas[d].total_counts = mcas[d].data.sum()
        else:
            print 'Unknown tag = '+tag+' in file: ' + file + '.'
    ##########
    fp.close()

    # Build dictionary to return
    r = {}
    r['n_detectors'] = n_detectors
    r['mcas'] = mcas
    r['rois'] = rois
    r['max_rois'] = max_rois
    r['environment'] = environment
    return r

#########################################################################
def write_file(detector, file):
    """
    Writes Mca or Med objects to a disk file.
    
    Inputs:
        detector: An med or mca object 
        file: The name of the disk file to write.
            
    Example:
        write_file(mca,'mca.001')
    """
    # Make sure detector is type MED for simplicity...
    try:
        if detector.det_type != "MED":
            if has_attr(detector,'environment'):
                env = detector.environment
                del detector.environment
                #delattr(detector,'environment')
            else:
                env = None
            detector = Med.Med(mcas=[detector])
            if env: detector.environment = env
    except:
        return
    
    write_ascii_file(detector, file)

#######################################################################
def write_ascii_file(med, file):
    """
    Writes Med data to a disk file.  The file 
    format is a tagged ASCII format.  The file contains the information 
    which it makes sense to store permanently, but 
    does not contain all of the internal state information for the detector.  
    Files written with this routine can be read with read_ascii_file().


    Inputs:
        med: Instance of Med object
        file: The name of the disk file to write.
    """
    # Get raw (uncorrected) data as a list, 
    data = []
    for mca in med.mcas:
        data.append(mca.data)
        if hasattr(mca,'rois'):
            rois.append(mca.rois)
    
    # Note we assume all mca data are the same length!
    # Also assume that mca.channels = [0.....len(mca.data)]
    # ie we dont write channels to file
    nchans = len(data[0])

    # Write header stuff
    n_det      = med.n_detectors
    start_time = med.mcas[0].start_time
    
    fformat = '%f ' * n_det
    eformat = '%e ' * n_det
    iformat = '%d ' * n_det
    sformat = '%s ' * n_det
    fp = open(file, 'w')
    fp.write('VERSION:    '+'3.1'+'\n')
    fp.write('ELEMENTS:   '+str(n_det)+'\n')
    fp.write('DATE:       '+str(start_time)+'\n')
    fp.write('CHANNELS:   '+str(nchans)+'\n')

    # count times and related    
    real_time=[]; live_time=[] ; input_counts=[]; tau=[]
    for mca in med.mcas:
        real_time.append(mca.real_time)
        live_time.append(mca.live_time)
        input_counts.append(mca.input_counts)
        tau.append(mca.tau)
    fp.write('REAL_TIME:    '+(fformat % tuple(real_time))+'\n')
    fp.write('LIVE_TIME:    '+(fformat % tuple(live_time))+'\n')
    fp.write('INPUT_COUNTS: '+(fformat % tuple(input_counts))+'\n')
    fp.write('TAU:          '+(fformat % tuple(tau))+'\n')

    # calibration data    
    offset=[]; slope=[]; quad=[]; two_theta=[]
    for mca in med.mcas:
        offset.append(mca.offset)
        slope.append(mca.slope)
        quad.append(mca.quad)
        two_theta.append(mca.two_theta)
    fp.write('CAL_OFFSET: '+(eformat % tuple(offset))+'\n')
    fp.write('CAL_SLOPE: '+(eformat % tuple(slope))+'\n')
    fp.write('CAL_QUAD: '+(eformat % tuple(quad))+'\n')
    fp.write('TWO_THETA: '+(fformat % tuple(two_theta))+'\n')

    # Write ROIS
    # note ROIS should always be in channel units!
    # Write number of rois for each mca
    nrois = []
    for d in range(n_det):
        if hasattr(med.mcas[d],'rois'):
            nrois.append(len(med.mcas[d].rois))
        else:
            nrois.append(0)
    fp.write('ROIS:       '+(iformat % tuple(nrois))+'\n')
    if max(nrois) > 0:
        for i in range(max(nrois)):
            num = str(i)
            left=[]; right=[]; label=[]
            for d in range(n_det):
                if (i < nrois[d]):
                    left.append(med.mcas[d].rois[i].left)
                    right.append(med.mcas[d].rois[i].right)
                    label.append(med.mcas[d].rois[i].label + '&')
                else:
                    left.append(0)
                    right.append(0)
                    label.append(' &')
            fp.write('ROI_'+num+'_LEFT:   '+(iformat % tuple(left))+'\n')
            fp.write('ROI_'+num+'_RIGHT:  '+(iformat % tuple(right))+'\n')
            fp.write('ROI_'+num+'_LABEL:  '+(sformat % tuple(label))+'\n')
    # Write environment
    if has_attr(med,'environment'):
        for e in med.environment:
            fp.write('ENVIRONMENT: '       + str(e.name) +
                                            '="'  + str(e.value) +
                                            '" (' + str(e.description) + ')\n')

    # Write data
    fp.write('DATA: \n')
    counts = Num.zeros(n_det)
    for i in range(nchans):
        for d in range(n_det):
            counts[d]=data[d][i]
        fp.write((iformat % tuple(counts))+'\n')

    # All done
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

    ext  = old_file[dot+1:]
    file = old_file[0:dot+1]
    nc   = str(len(ext))
    try:
        # Convert to number, add one, catch error
        ext = int(ext)+1      
        format = '%' + nc + '.' + nc + 'd'
        ext = (format % ext)
        new_file = file + ext
        return new_file
    except:
        return old_file

#########################################################################
#########################################################################
def test():
    return read_med('test.xrf')

if __name__ == "__main__":
    med = test()
    print med
    