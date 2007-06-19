# T. Trainor (2006)
# Wrapper functions for spec files
# T. Trainor, 6-10-2006
#
# --------------
# Modifications
# --------------
#
#
# --------------
# Todo
# --------------
import Util
from Num import Num
import os
import types
from SpecFile import SpecFile
#import SD as ScanData
import ScanData
import XRF

####################################################
def read_spec(fname,tdl=None,**kws):
    """    >>sf = read("file_name")
    >>spec.read file_name
    
    The spec file variable is an object that contains all the information
    this can be passed to the various functions for extracting and plotting
    scan data. Note that the spec file_name can contain a path.  The path 
    may be full, relative to the current directory or relative to _sys.work.
    When using command syntax the variable name will be the file name prefix
    and stored in the spec group.
    """
    sf = SpecFile(fname)
    return sf

def read_spec_cmd(val,tdl=None,**kws):
    name = val._fname
    if '.' in name: name = name.split('.',1)[0]
    name = 'spec.%s' % name
    tdl.setVariable(name,val=val)

def show_scan(sf,scan=None,all=False,**kws):
    """    >>spec.show spec_file, [scan=10,all=True]
    
    Show spec file/scan information.  If no scan number is provided this outputs
    a summary of the data file other wise a summary of the scan is provided.
    The all flag controls the level of detail.  The spec_file argument may be a file name
    of a spec_file object (see read_spec).  This function has no return value. """
    
    if type(sf) == types.StringType:
        sf = SpecFile(sf)
    if not sf._ok: return
    
    if scan == None:
        print sf
        return
    if sf._check_range(scan):
        sum = sf.scan_info(scan)
    else:
        sum = None
    if sum:
        #keys = sum.keys()
        #keys.sort()
        #for k in keys:
        #    print "#%s  %s\n" % (k, sum[k])
        print "#Scan:     %i" % sum['index']
        print "#Date:     %s" % sum['date'].strip()
        print "#Command:  %s" % sum['cmd'].strip()
        print "#ATTEN:    %s" % sum['atten'].strip()
        print "#Labels:   %s" % sum['labels'].strip()
        if all:
            print "#G:"
            print Util.show_list(sum['G'].split())
            print "#P:"
            m_names = sum['mot_names'].split()
            p_vals  = sum['P'].split()
            if len(m_names) != len(p_vals):
                print "Mismatch in Motor Names and Motor Values"
            else:
                lout = []
                for j in range(len(m_names)):
                    x = "%s=%s" % (m_names[j], p_vals[j])
                    lout.append(x)
                print Util.show_list(lout)        
    else:
        print "Scan not found"
        print sf

def scan_data(sf,scan,**kws):
    """    >>dat = spec.data(spec_file,scan)
    
    The spec_file argument may be a file name or a spec_file object (see read_spec).
    The returned array is the scan data:
    dat[0] = 1st column (as a row vector)
    dat[1] = 2nd column (as a row vector) etc..
    """
    if type(sf) == types.StringType:
        sf = SpecFile(sf)
    if not sf._ok: return None
    d = Num.array(sf.scan_data(scan))
    return d.transpose()

def scan_cols(sf,scan,cols=None,**kws):
    """    >>dat = spec.col(spec_file,scan,cols=[column_lables])
    
    The spec_file argument may be a file name or a spec_file object (see read_spec).
    The returned array is the scan data:
    dat[0] = column mathing the first column label (or the first column - default)
    dat[1] = column mathing the second column label (or the last column - default)
    """
    if type(sf) == types.StringType:
        sf = SpecFile(sf)
    if not sf._ok: return None
    d = sf.scan_dict(scan)
    dd = d['data']
    lbls = d['labels']
    ncol = d['ncol']
    dat = []
    if cols:
        if type(cols) != types.ListType: cols = [cols]
        for c in cols:
            if c in dd.keys():
                dat.append(dd[c])
            else:
                print "Warning: skipping column label '%s', not found" % c
    else:
        dat.append(dd[lbls[0]])
        dat.append(dd[lbls[ncol-1]])
    return dat


def scan_dict(sf,scan,**kws):
    """    >>dat = spec.dict(spec_file,scan)
    
    The spec_file argument may be a file name or a spec_file object (see read_spec).
    The returned dictionary contains all the scan data and relevant scan information.
    Use >dictkeys(dat) to veiw the contents.
    """
    if type(sf) == types.StringType:
        sf = SpecFile(sf)
    if not sf._ok: return None
    return sf.scan_dict(scan)

def get_scan(sf,scan):
    """    >>dat = spec.get_scan(spec_file,scan)
    
    The spec_file argument may be a file name or a spec_file object (see read_spec).
    The returned value is a ScanData class.
    """
    if type(sf) == types.StringType:
        sf = SpecFile(sf)
    if not sf._ok: return None
    return sf.get_scan(scan)

def filename(sf,**kws):
    """    >>file = spec.filename(spec_file)
    
    return the filename
    """
    if not sf._ok: return None
    return sf._fname


#################################################
_groups_ = [('spec',True)]
_func_ = {"spec.read":(read_spec,read_spec_cmd),
          "spec.show":(show_scan,None),
          "spec.data":scan_data,
          "spec.col":scan_cols,
          "spec.dict":scan_dict,
          "spec.get_scan":get_scan,
          "spec.filename":filename}

_scripts_ = ['spec_scripts.tdl']
