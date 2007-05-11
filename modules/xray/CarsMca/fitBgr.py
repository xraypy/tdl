## Automatically adapted for numpy.oldnumeric May 08, 2007 by 

########################################################################
# From Mca
########################################################################

import copy
import numpy.oldnumeric as Numeric
import CARSMath

########################################################################
class McaBackground:
    """
    Defines the parameters for fitting backgrounds in Mca objects.
    Fields and default values:
        .slope        = 1.0
        .exponent     = 2
        .top_width    = 0.
        .bottom_width = 4.
        .tangent      = 0
        .compress     = 4
        
    See the documentation on fit_background() for information on the meaning
    of these fields.
    """
    def __init__(self,slope=1.0, exponent=2,top_width=0.0,bottom_width=4.0,tangent=0,compress=4):
        self.slope        = slope
        self.exponent     = exponent
        self.top_width    = top_width
        self.bottom_width = bottom_width
        self.tangent      = tangent
        self.compress     = compress


def fit_background(data, params):
    """
    This function fits a background to an MCA spectrum. The background is
    fitted using an enhanced version of the algorithm published by
    Kajfosz, J. and Kwiatek, W .M. (1987)  "Non-polynomial approximation of
    background in x-ray spectra." Nucl. Instrum. Methods B22, 78-81.

    Inputs:
        data: an array of count values, eg from an mca spectrum.  
    
        params: (should hold the following)
             slope:
                Slope for the conversion from channel number to energy.
                Ie the slope from calibration
                
             top_width:
                Specifies the width of the polynomials which are concave upward.
                The top_width is the full width in energy units at which the
                magnitude of the polynomial is 100 counts. The default is 0, which
                means that concave upward polynomials are not used.

             bottom_width:
                Specifies the width of the polynomials which are concave downward.
                The bottom_width is the full width in energy units at which the
                magnitude of the polynomial is 100 counts. The default is 4.

             exponent:
                Specifies the power of polynomial which is used. The power must be
                an integer. The default is 2, i.e. parabolas. Higher exponents,
                for example EXPONENT=4, results in polynomials with flatter tops
                and steeper sides, which can better fit spectra with steeply
                sloping backgrounds.

             tangent:
                Specifies that the polynomials are to be tangent to the slope of the
                spectrum. The default is vertical polynomials. This option works
                best on steeply sloping spectra. It has trouble in spectra with
                big peaks because the polynomials are very tilted up inside the
                peaks.

             compress:
                Compression factor to apply before fitting the background.
                Default=4, which means, for example, that a 2048 channel spectrum
                will be rebinned to 512 channels before fitting.
                The compression is done on a temporary copy of the input spectrum,
                so the input spectrum itself is unchanged.
                The algorithm works best if the spectrum is compressed before it
                is fitted. There are two reasons for this. First, the background
                is constrained to never be larger than the data itself. If the
                spectrum has negative noise spikes they will cause the fit to be
                too low. Compression will smooth out such noise spikes.
                Second, the algorithm requires about 3*N^2 operations, so the time
                required grows rapidly with the size of the input spectrum. On a
                200 MHz Pentium it takes about 3 seconds to fit a 2048 channel
                spectrum with COMPRESS=1 (no compression), but only 0.2 seconds
                with COMPRESS=4 (the default).
                
 Procedure:
         1) At each channel "i" an n'th degree polynomial which is concave up
         is fitted. Its equation is

                                      n
                         (e(i) - e(j))
         f(j,i) = y(i) + --------------
                                       n
                            top_width

         where f(j,i) is the fitted counts in channel j for the polynomial
         centered in channel i. y(i) is the input counts in channel "i", e(i) is
         the energy of channel i, e(j) is the energy of channel j, and
         "top_width" and "n" are user-specified parameters. The background count
         in channel "j", b(j) is defined as

         b(j) = min ((f(j,i), y(j))
                 i

         b(j) is thus the smallest fitted polynomial in channel j, or the raw
         data, whichever is smaller.

         2) After the concave up polynomials have been fitted, a series of
         concave down polynomials are constructed. At each channel "i" an n'th
         degree polynomial which is concave up is fitted. The polynomial is slid
         up from below until it "just touches" some channel of the spectrum. Call
         this channel "i". The maximum height of the polynomial is thus

                                                n
                                    (e(i) - e(j))
         height(j) = max ( b(j) +  --------------  )
                        i                        n
                                    bottom_width

         where bottom_width is a user_specified parameter.

         3) Once the value of height(i) is known the polynomial is fitted. The
         background counts in each channel are then determined from:

                                                 n
                                    (e(i) - e(j))
         bgd(j) = max ( height(i) + --------------
                    i                             n
                                     bottom_width

         bgd(j) is thus the maximum counts for any of the concave down
         polynomials passing though channel j.

         Before the concave-down polynomials are fitted the spectrum at each
         channel it is possible to subtract out a straight line which is
         tangent to the spectrum at that channel. Use the /TANGENT qualifier to
         do this. This is equivalent to fitting a "tilted" polynomial whose
         apex is tangent to the spectrum at that channel. By fitting
         polynomials which are tangent rather than vertical the background fit
         is much improved on spectra with steep slopes.

    Outputs:
         This function returns the optimized background for data.

    Example:
        bgd = fit_background(mca, params)
    """
    
    REFERENCE_AMPL=100.
    TINY = 1.E-20
    HUGE = 1.E20
    MAX_TANGENT=2

    slope        = params.slope
    exponent     = params.exponent
    top_width    = params.top_width
    bottom_width = params.bottom_width
    tangent      = params.tangent
    compress     = params.compress

    bgd         = copy.copy(data)
    nchans      = len(data)
    scratch     = copy.copy(bgd)

    # Compress scratch spectrum
    if (compress > 1):
        scratch = CARSMath.compress_array(scratch, compress)
        slope = slope * compress
        nchans = nchans / compress

    # Copy scratch spectrum to background spectrum
    bckgnd = copy.copy(scratch)

    # Find maximum counts in input spectrum. This information is used to
    # limit the size of the function lookup table
    max_counts = max(scratch)

    #  Fit functions which come down from top
    if (top_width > 0.):
        #   First make a lookup table of this function
        chan_width = top_width / (2. * slope)
        denom = chan_width**exponent
        indices = Numeric.arange(float(nchans*2+1)) - nchans
        power_funct = indices**exponent * (REFERENCE_AMPL / denom)
        power_funct = Numeric.compress((power_funct <= max_counts), power_funct)
        max_index = len(power_funct)/2 - 1

        for center_chan in range(nchans):
            first_chan = max((center_chan - max_index), 0)
            last_chan = min((center_chan + max_index), (nchans-1))
            f = first_chan - center_chan + max_index
            l = last_chan - center_chan + max_index
            test = scratch[center_chan] + power_funct[f:l+1]
            sub = bckgnd[first_chan:last_chan+1] 
            bckgnd[first_chan:last_chan+1] = Numeric.maximum(sub, test)

    # Copy this approximation of background to scratch
    scratch = copy.copy(bckgnd)

    # Find maximum counts in scratch spectrum. This information is used to
    #   limit the size of the function lookup table
    max_counts = max(scratch)

    # Fit functions which come up from below
    bckgnd = Numeric.arange(float(nchans)) - HUGE

    # First make a lookup table of this function
    chan_width = bottom_width / (2. * slope)
    if (chan_width == 0.):
        denom = TINY
    else:
        denom = chan_width**exponent
    
    indices = Numeric.arange(float(nchans*2+1)) - nchans
    power_funct = indices**exponent  * (REFERENCE_AMPL / denom)
    power_funct = Numeric.compress((power_funct <= max_counts), power_funct)
    max_index = len(power_funct)/2 - 1

    for center_chan in range(nchans-1):
        tangent_slope = 0.
        if (tangent):
            # Find slope of tangent to spectrum at this channel
            first_chan = max((center_chan - MAX_TANGENT), 0)
            last_chan = min((center_chan + MAX_TANGENT), (nchans-1))
            denom = center_chan - Numeric.arange(float(last_chan - first_chan + 1))
            tangent_slope = (scratch[center_chan] - 
                             scratch[first_chan:last_chan+1]) / max(denom, 1)
            tangent_slope = Numeric.sum(tangent_slope) / (last_chan - first_chan)

        first_chan = max((center_chan - max_index), 0)
        last_chan = min((center_chan + max_index), (nchans-1))
        last_chan = max(last_chan, first_chan)
        nc = last_chan - first_chan + 1
        lin_offset = scratch[center_chan] + \
                     (Numeric.arange(float(nc)) - nc/2) * tangent_slope

        # Find the maximum height of a function centered on this channel
        # such that it is never higher than the counts in any channel

        f = first_chan - center_chan + max_index
        l = last_chan - center_chan + max_index
        test = scratch[first_chan:last_chan+1] - lin_offset + \
                                                 power_funct[f:l+1]
        height = min(test)

        # We now have the function height. Set the background to the
        # height of the maximum function amplitude at each channel

        test = height + lin_offset - power_funct[f:l+1]
        sub = bckgnd[first_chan:last_chan+1]
        bckgnd[first_chan:last_chan+1] = Numeric.maximum(sub, test)

    # Expand spectrum
    if (compress > 1):
        bckgnd = CARSMath.expand_array(bckgnd, compress)
    bgd = bckgnd.astype(Numeric.Int)

    return bgd

