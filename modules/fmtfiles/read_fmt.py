"""
Read simple formatted files

Authors/Modifications:
----------------------
Tom Trainor (trainor@alaska.edu) 


"""
##########################################################################

import numpy as num

##########################################################################
def read_column(fname): 
    """
    Read column data from file
    
    Parameters:
    -----------
    * fname is the filename to read

    Returns:
    --------
    * array of data from file

    Examples:
    ---------
    >> data = read_column('data.dat')

    data.dat contents:
    # comments
    1 2 3
    4 5 6
    7 8 9

    """
    f = open(fname)
    lines = f.readlines()
    f.close()
    data = []
    for line in lines:
        if line[0] == '#':
            pass
        else:
            tmp = line.split()
            #print tmp
            d = map(float,tmp)
            #print d
            data.append(d)
    data = num.array(data)
    return num.transpose(data)

