'''
Filter Tools
Author: Craig Biwer (cbiwer@uchicago.edu)
2/13/2012
'''

import h5py
import numpy

# Given two dates, checks that the first comes before the second
# Format: 'Day Month Date Time Year', e.g. 'Mon Jan 12 08:42:32 2011'
def is_before(am_i, before_me):
    all_months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
                  'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
    (test_day, test_month, test_date, test_time, test_year) = am_i.split()
    (other_day, other_month, other_date,
                             other_time, other_year) = before_me.split()
    if int(test_year) < int(other_year):
        return True
    elif int(test_year) == int(other_year):
        if all_months.index(test_month) < all_months.index(other_month):
            return True
        elif all_months.index(test_month) == all_months.index(other_month):
            if int(test_date) < int(other_date):
                return True
            elif int(test_date) == int(other_date):
                if test_time < other_time:
                    return True
    return False
    
# Given two dates, checks that the first comes after the second
# Format: 'Day Month Date Time Year', e.g. 'Mon Jan 12 08:42:32 2011'
def is_after(am_i, after_me):
    return am_i == after_me or not is_before(am_i, after_me)

# Given a series of lists, returns the union of those lists
def list_union(*args):
    for entry in args:
        if entry is None:
            args = list(args)
            args.remove(None)
            return list_union(*args)
    union = [item for sublist in args for item in sublist]
    union = list(set(union))
    return union

# Given a series of lists, returns the intersection of those lists
def list_intersect(*args):
    intersect = set(args[0])
    for entry in args:
        intersect = intersect.intersection(entry)
    return list(intersect)

# The heart of the filtering process
# Given an opened h5 file, filters on the given criteria
# and returns a list of the matched scans / points
def cases(in_here, of_this, such_that, this_level='scan'):
    '''
    Filters the given file or list of scans and returns those that pass the
    specified criteria. Returns None if an error is encountered.

    in_here is the opened h5 file [f = h5py.File(<filename>, 'r')]
    of_this is the attribute / label to filter
    such_that is the criteria you want the filtered variable to pass
        eg '< 1.4' or '== "rodscan"' (note the quotes: this is passed to eval)
    this_level is the sensitivity of the filter:
        scan: default; checks the attribute / label of_this against such_that
        point: checks all the points and returns only those that match

    If this_level is scan, returns a list of paths in the h5 file
    If this_level is point, returns a list of (path, index) tuples
    '''
    
    if isinstance(in_here, h5py.File):
        all_possible = []
        for spec, group in in_here.items():
            for number, scan in group.items():
                all_possible.append(scan.name)
    else:
        print 'Error: input file not recognized'
        return None
    
    return_this = []
    if this_level == 'scan':
        for scan in all_possible:
            this_ind = None
            this_value = None
            if in_here[scan].attrs.get(of_this, None) is not None:
                this_value = in_here[scan].attrs.get(of_this)
            elif of_this in in_here[scan]['param_labs']:
                this_ind = list(in_here[scan]['param_labs']).index(of_this)
                this_value = in_here[scan]['param_data'][this_ind]
            '''
            else:
                print 'Error: Unrecognized attribute or label'
                return None
            '''
            if this_value is None:
                pass
            elif isinstance(this_value, basestring):
                if eval('"' + this_value + '" ' + such_that):
                    return_this.append(scan)
            elif isinstance(this_value, (int, float, bool, numpy.bool_)):
                if eval(str(this_value) + " " + such_that):
                    return_this.append(scan)
            else:
                print 'Error: Unrecognized type'
                print this_value
                return None
        return return_this
                
    elif this_level == 'point':
        for scan in all_possible:
            this_ind = None
            this_value = None
            if in_here[scan].attrs.get(of_this, None) is not None:
                this_value = in_here[scan].attrs.get(of_this)
                if isinstance(this_value, basestring):
                    if eval('"' + this_value + '" ' + such_that):
                        for i in range(len(in_here[scan]['point_data'])):
                            return_this.append((scan, i))
                elif isinstance(this_value, (int, float, bool, numpy.bool_)):
                    if eval(str(this_value) + " " + such_that):
                        for i in range(len(in_here[scan]['point_data'])):
                            return_this.append((scan, i))
                else:
                    print 'Error: Unrecognized type'
                    print this_value
                    return None
            elif of_this in in_here[scan]['point_labs']:
                this_ind = list(in_here[scan]['point_labs']).index(of_this)
                for i in range(len(in_here[scan]['point_data'])):
                    this_value = in_here[scan]['point_data'][i][this_ind]
                    if isinstance(this_value, basestring):
                        if eval('"' + this_value + '" ' + such_that):
                            return_this.append((scan, i))
                    elif isinstance(this_value,
                                    (int, float, bool, numpy.bool_)):
                        if eval(str(this_value) + " " + such_that):
                            return_this.append((scan, i))
                    else:
                        print 'Error: Unrecognized type'
                        print this_value
                        return None
            elif of_this in in_here[scan]['param_labs']:
                this_ind = list(in_here[scan]['param_labs']).index(of_this)
                this_value = in_here[scan]['param_data'][this_ind]
                if isinstance(this_value, basestring):
                    if eval('"' + this_value + '" ' + such_that):
                        for i in range(len(in_here[scan]['point_data'])):
                            return_this.append((scan, i))
                elif isinstance(this_value, (int, float, bool, numpy.bool_)):
                    if eval(str(this_value) + " " + such_that):
                        for i in range(len(in_here[scan]['point_data'])):
                            return_this.append((scan, i))
                else:
                    print 'Error: Unrecognized type'
                    print this_value
                    return None
        return return_this
    else:
        print 'Error: unrecognized filter level'
        return None