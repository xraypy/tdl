########################################################################
"""
Mark Rivers, GSECARS
MultiChannel Analyzer (MCA) Energy calibration.

Modifications:
--------------
- See http://cars9.uchicago.edu/software/python/index.html
- Modified for Tdl, tpt

"""
########################################################################
"""
Compute channel to energy conversions and visa-versa
for multi-channel analyzer.

We use the following convention:

    energy = offset + channels * slope

where channels is an integer array (energy is a double array).
Note that assumed that the offset is ALWAYS defined to be
relative to the channel number of zero

Therefore, the index of the channel array does not
need to equal its value.  E.g. 

  channels = num.arange(nchan,dtype=int)        ==> [0,1,2................2048]
  channels = num.arange(mi_ch,ma_ch,dtype=int)  ==> [33,34,..............657]

are both valid.  

Note make shorter names for the functions:
 - ch2en = channel_to_energy
 - en2ch = energy_to_channel

"""
########################################################################

import numpy as num

########################################################################
def channel_to_energy(channels, offset=0.0, slope=1.0, quad=0.0):
    """
    Converts channels to energy using the current calibration values
    for the Mca.  

    Inputs:
        channels: The channel numbers to be converted to energy.  This can be
                  a single number or a sequence of channel numbers.
            
    Outputs:
        This function returns the equivalent energy for the input channels.
    """
    c  = num.asarray(channels,dtype=num.double)
    if quad != 0.:
        return offset +  c * (slope + c * quad)
    else:
        return offset +  c * (slope)

#########################################################################
def energy_idx(energy,emin=-1.,emax=-1.):
    """
    Get the indicies of the energy array for emin and emax:

    idx = energy_idx(3.2,7.1,energy)
    en  = energy[idx]
    dat = data[idx]
    """
    # find emin
    if emin >= 0.0:
        delta    = num.abs(energy - emin)
        idx_min  = num.where( delta == min(delta) )
        idx_min  = idx_min[0][0]
    else:
        idx_min = 0

    # find emax
    if emax >= 0.0 and emax >= emin:
        delta    = num.abs(energy - emax)
        idx_max  = num.where( delta == min(delta) )
        idx_max  = idx_max[0][0]
    else:
        idx_max = len(energy) - 1
 
    #return (idx_mi,idx_ma)
    index = num.arange(idx_min,idx_max+1,dtype=num.int)
    return index

########################################################################
def energy_to_channel(energy, offset=0.0, slope=1.0, quad=0.0, clip=0):
    """
    Converts energy to channel numbers for an Mca.  

    Inputs:
        energy: The energy values to be converted to channels. This can be a
                single number or a sequence energy values.
            
    Keywords:
        clip: Set this flag to >0 to clip the returned values to be between
              0 and clip (inclusive).  The default is not to clip.
            
    Outputs:
        This function returns the closest equivalent channel for the
        input energy.  
        
    """
    if (quad == 0.0):
        channel = ((energy - offset) / slope)
    else:
        # Use the quadratic formula
        a = quad
        b = slope
        c = offset - energy
        # There are 2 roots.  I think we always want the "+" root?
        channel = (-b + num.sqrt(b**2 - 4.*a*c))/(2.*a)
    channel = num.around(channel)

    if (clip > 0):
        channel = num.clip(channel, 0, clip)
        # Note if energy is an array, below may
        # result in a return array of different length
        #condition = (channel >= 0) and (channel <=clip)
        #channel = channel.compress(condition)

    channel = channel.astype(num.int)

    return channel

########################################################################
def channel_to_d(channels,two_theta=0.0, offset=0.0, slope=1.0, quad=0.0):
    """
    Converts channels to "d-spacing" 

    Inputs:
        channels: The channel numbers to be converted to "d-spacing".
                  This can be a single number or a list of channel numbers.
            
    Outputs:
        This function returns the equivalent "d-spacing" for the input channels.
        The output units are in Angstroms.
        
    Restrictions:
        This function assumes that the units of the energy calibration are keV
        and that the units of "two-theta" are degrees.
        
    Example:
        channels = [100,200,300]
        d = mca.channel_to_d(channels)
    """
    e = channel_to_energy(channels,offset=offset, slope=slope, quad=quad)
    return 12.398 / (2. * e * num.sin(two_theta*num.pi/180.))


########################################################################
def d_to_channel(d, two_theta=0.0, offset=0.0, slope=1.0, quad=0.0, clip=0):
    """
    Converts "d-spacing" to channels

    Inputs:
        d:  The "d-spacing" values to be converted to channels.
            This can be a single number or an array of values.
            
    Keywords:
        clip: Set this flag to 1 to clip the returned values to be between
              0 and nchans-1.  The default is not to clip.
            
    Outputs:
        This function returns the closest equivalent channel for the input
        "d-spacing". 
        
    Example:
        channel = d_to_channel(1.598)
    """
    e = 12.398 / (2. * d * num.sin(two_theta*num.pi/180./2.))
    return energy_to_channel(e,offset=offset, slope=slope, quad=quad, clip=clip)

########################################################################
########################################################################
if __name__ == "__main__":
    e = num.arange(1,15.,.1)
    idx = energy_idx(e,1.66,1.66)
    print e[idx]
    
