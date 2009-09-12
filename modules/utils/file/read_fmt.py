##########################################################################
"""
Tom Trainor (fftpt@uaf.edu ) 

Modifications:
--------------


"""
##########################################################################

import numpy as num

##########################################################################
def read_col_data(fname): 
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
