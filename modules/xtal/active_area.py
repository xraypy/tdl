"""
Active area calculations

Authors/Modifications:
----------------------
Tom Trainor (tptrainor@alaska.edu)


Notes:
------
These area are based on completley geometric projections.  ie we assume
that the beam is perfectly parrallel (no divergence)

Note would be interesting if we could set the detector polygon from an 
arbirary set of pixels (e.g. roi) on image detector... especially if 
shifted off center (or better off center relative to a specific hkl)
"""
##########################################################################

import numpy as num
import types
from matplotlib import pyplot

from mathutil import cosd, sind, tand
from mathutil import arccosd, arcsind, arctand
from mathutil import cartesian_mag, cartesian_angle

from polygon import poly_area, poly_area_num, inner_polygon
from polygon import plot_polygon, plot_points, plot_circle

##########################################################################
def active_area(nm,ki=num.array([0.,1.,0.]),kr=num.array([0.,1.,0.]),
                beam=[],det=None,sample=1.,plot=False,fig=None):
    """
    Calc the area of overlap of beam, sample and detector
    surface polygon projections.  Return:
        A_beam = total area of the beam projection into the surface plane
        A_int = area of sample illuminated within the detector projection 
    Use to correct scattering data for area effects, including spilloff, i.e.
        A_ratio = A_int/A_beam 
        Ic = I/A_ratio,   I = Idet/Io

    nm is the surface normal vector, defined in the laboratory frame.

    ki and kr are vectors defined in the lab frame parallel to the
    incident beam and diffracted beam directions respectively
    (magnitudes are arbitrary)

    beam, det and sample are lists that hold the the lab-frame (3D) 
    vectors defining the beam apperature, detector apperature
    and sample shape 
    
    If det = None then we ignore it and just compute spill-off correction

    If sample is a single number we take it as the diamter of a round sample
    mounted flat. Otherwise sample hould be a list of lab frame vectors that 
    describes the sample polygon.

    If sample == None, then we assume the sample is infinite in size    
    
    The lab frame coordinate systems is defined such that:
        x is vertical (perpendicular, pointing to the ceiling of the hutch)
        y is directed along the incident beam path
        z make the system right handed and lies in the horizontal scattering plane
          (i.e. z is parallel to the phi axis)
    The center (0,0,0) of the lab frame is the rotation center of the instrument.

    If plot = True then makes plot

    """
    # Calc surf system transformation matrix
    M = calc_surf_transform(nm)

    # calc k vectors in the surf frame (i.e. [xs,ys,zs])
    # first make them unit vectors, even though no real need to..
    ki = ki/cartesian_mag(ki)
    kr = kr/cartesian_mag(kr)
    ki_s = num.dot(M,ki)
    kr_s = num.dot(M,kr)

    #################################################################
    # Get the beam, detector and surface polygons.
    # Convert vectors to surface frame.
    # Then for beam and det project the vectors onto
    # the surface along the k vectors.  ie. calc the vector
    # intercepts with the surface plane
    # Each of the below polygons is a list of 2D vectors
    # ie in-plane [x_s, y_s] points 
    #################################################################
    # beam 
    if type(beam) == types.ListType:
        beam_poly = []
        if len(beam) < 3:
            print "Error in beam description"
            return None
        for v in beam:
            vs  = num.dot(M,v)
            int = surface_intercept(ki_s, vs)
            beam_poly.append(int)
    else:
        beam_poly=None
    # det 
    if type(det) == types.ListType:
        det_poly  = []
        if len(det) < 3:
            print "Error in det description"
            return None
        for v in det:
            vs  = num.dot(M,v)
            int = surface_intercept(kr_s, vs)
            det_poly.append(int)
    else:
        det_poly=None
    # sample 
    if type(sample) == types.ListType:
        sam_poly  = []
        if len(sample) < 3:
            print "Error in sample description"
            return None
        for v in sample:
            vs  = num.dot(M,v)
            if num.fabs(vs[2]) > 0.01:
                print "Warning sample hieght problem"
            sam_poly.append(vs[:2])
        sample_shape=True
    else:
        sam_poly = None
        sample_shape = False

    #####################################################################
    # depending on whether we have a sample shape description
    # or assume a round sample, pass along appropriate
    # data to compute areas.  Note all the polygons are lists
    # of 2D surface frame vectors (ie in plane vectors)
    #####################################################################
    if sample_shape == False:
        (A_beam,A_int) = _area_round(beam_poly,det_poly,diameter=sample,plot=plot,fig=fig)
    else:
        (A_beam,A_int) = _area_polygon(beam_poly,det_poly,sam_poly,plot=plot,fig=fig)
    return (A_beam,A_int)

#########################################################################
def _area_round(beam_poly,det_poly,diameter=None,plot=False,fig=None):
    """
    compute areas for round sample of fixed diameter
    if det_poly = None, then just compute the beam and
    sample overlap ie A_int/A_beam = spill fraction
    """
    if plot:
        if fig != None: pyplot.figure(fig)
        pyplot.clf()
    A_beam = poly_area(beam_poly)
    if det_poly != None:
        inner_poly = inner_polygon(beam_poly,det_poly)
    else:
        inner_poly = beam_poly
    if diameter <= 0.:
        diameter = None
    if diameter == None:
        A_int = poly_area(inner_poly)
    else:
        A_int = poly_area_num(beam_poly,diameter=diameter,plot=plot)

    # make plots
    if plot:
        plot_polygon(beam_poly,fmt='ro-',label='Beam')
        if det_poly != None:
            plot_polygon(det_poly,fmt='ko-',label='Detector')
            plot_points(inner_poly,fmt='go-')
            plot_polygon(inner_poly,fmt='g--',linewidth=4,label='Intersection')
        if diameter!= None:
            plot_circle(diameter/2.,fmt='b-',label='Sample')
        pyplot.xlabel('Xs')
        pyplot.ylabel('Ys')
        ll = round(num.sqrt(4.*A_int))
        pyplot.xlim(-ll,ll)
        pyplot.ylim(-ll,ll)
        pyplot.grid()
        pyplot.legend()

    return (A_beam,A_int)

#########################################################################
def _area_polygon(beam_poly,det_poly,sam_poly,plot=False,fig=None):
    """
    compute areas for polgon sample
    if det_poly = None, then just compute the beam and
    sample overlap ie A_int/A_beam = spill fraction
    """
    if plot:
        if fig != None: pyplot.figure(fig)
        pyplot.clf()
    A_beam = poly_area(beam_poly)
    if sam_poly == None:
        inner_poly = beam_poly
    else:
        inner_poly = inner_polygon(beam_poly,sam_poly)
    if det_poly != None:
        inner_poly = inner_polygon(inner_poly,det_poly)
    A_int = poly_area(inner_poly)
    
    # make plots
    if plot:
        plot_polygon(beam_poly,fmt='ro-',label='Beam')
        if sam_poly != None:
            plot_polygon(sam_poly,fmt='bo-',label='Sample')
        if det_poly != None:
            plot_polygon(det_poly,fmt='ko-',label='Detector')
        plot_points(inner_poly,fmt='go-')
        plot_polygon(inner_poly,fmt='g--',linewidth=4,label='Intersection')
        pyplot.xlabel('Xs')
        pyplot.ylabel('Ys')
        ll = round(num.sqrt(4.*A_int))
        pyplot.xlim(-ll,ll)
        pyplot.ylim(-ll,ll)
        pyplot.grid()
        pyplot.legend()

    return (A_beam,A_int)

##########################################################################
def calc_surf_transform(nm):
    """
    This routine calculates the matrix, M, which transforms the
    indicies of vectors defined in the lab frame basis to the
    a surface frame defined such that:
        zs is along the surface normal (parrallel to nm)
        ys is the projection of the lab frame -y axis onto the surface
        xs is defined to make it a right handed orthonormal set
           (and is therefore in the surface plane)

    nm is the surface normal vector, defined in the 
    laboratory m-frame (i.e. pointing in an arbitrary
    direction for a particular gonio setting)
    - see gonio_psic.py for more details.

    The surface transform is defined as follows:
        |e_xs|      |e_x|              |vx|     
        |e_ys| = F* |e_y|    and   F = |vy|   
        |e_zs|      |e_z|              |vz|    
    where e's are the (cartesian) basis vectors.  Given a
    vector [x,y,z] defined in the lab frame basis [e_x,e_y,e_z],
    we can compute the indicies of this vector in the surface
    basis [e_xs,e_ys,e_zs] from:
        |xs|      |x|
        |ys| = M* |y|
        |zs|      |z|
    where
        M = transpose(inv(F)) = inv(transpose(F))
        
    Note see lattice.py for more notes on general transforms, and see
    gonio_psic.py for notes on calc of the surface normal vector.  
    """
    v_z = nm / cartesian_mag(nm)

    v   = num.array([0.,-1.,0.])
    v_y = v - (num.dot(v,v_z) * v_z) / (cartesian_mag(v_z)**2.) 
    v_y = v_y / cartesian_mag(v_y)
    
    v_x = num.cross(v_y,v_z)
    v_x = v_x / cartesian_mag(v_x)

    F = num.array([v_x,v_y,v_z])
    M = num.linalg.inv(F.transpose())

    return M

##################################################################
def surface_intercept(k,v):
    """
    Find surface coordinates [x,y] when the tip of vector v is
    projected into the surface plane along k.  ie [x,y] for z=0
    along k passing through the point defined by v
    """
    # if k[2] = 0 then k is parallel to xy plane
    # and no intercept is possible.  
    if k[2] == 0: return None
    
    x_int = v[0] - (k[0]/k[2])*v[2]
    y_int = v[1] - (k[1]/k[2])*v[2]

    vi = num.array([x_int, y_int])

    return vi

##########################################################################
def surface_intercept_bounds(k,v,diameter):
    """
    Find z-intercepts within bounds
    
    We assume that the plane is a circle of fixed diameter, 
    and scale back the intercept to the edge of the circle,
    This is done along by subtracting a vector that is parrelel
    to the in-plane projection of the k vector. The magnitude of
    the subtracted vector is set so the result has magnitude
    equal to the diameter.
    
    If that doesnt result in an intercept then just rescale
    the original intercept vector so its magnitude is the
    diameter.
    
    This is an approximate way to handle intercepts for
    circular samples, not recommended for use(!)
    """
    vi = surface_intercept(k,v)
    r = cartesian_mag(vi)
    if r <= diameter/2.:
        return vi
    # if r>diameter/2 need to scale it back
    # Therefore solve the following for mag
    #    vi_new = vi - ks*mag/norm(ks),
    #    mag(vi_new) = diameter/2
    # therefore end up with a quadratic:
    ks = num.array(k[0:2])
    a = 1.
    b = -2. * (v[0]*ks[0] + vi[1]*ks[1])/cartesian_mag(ks)
    c = cartesian_mag(vi)**2. - (diameter/2.)**2.
    arg = b**2. - 4.*a*c
    # make sure this method will result in a solution
    if arg > 0.: 
        mag1 = (-b + num.sqrt(arg))/(2.*a)
        mag2 = (-b - num.sqrt(arg))/(2.*a)
        vi1 = vi - ks*mag1/cartesian_mag(ks)
        vi2 = vi - ks*mag2/cartesian_mag(ks)
        # chose soln that makes smallest angle 
        # with the original intercept
        angle1 = num.fabs(cartesian_angle(vi,vi1))
        angle2 = num.fabs(cartesian_angle(vi,vi2))
        if angle1 < angle2:
            return vi1
        else:
            return vi2
    # otherwise just rescale the original vector
    else:
        vi = vi * (diameter/2.) / cartesian_mag(vi)
        return vi

##########################################################################
##########################################################################
def test1():
    r = 0.5
    k = [0., 1., -0.1]
    v = [.2, 0.,  1.]
    v1 = surface_intercept(k,v)
    v2 = surface_intercept_bounds(k,v,2.*r)
    print v1,v2
    plot_circle(r)
    plot_points([v1,v2])
    pyplot.grid()

def test2():
    import gonio_psic
    psic = gonio_psic.test2(show=False)
    psic.set_angles(phi=42.,chi=33,eta=20.,
                    mu=15.,nu=75.,delta=20.)
    # get beam and detector vectors
    beam = gonio_psic.beam_vectors(h=1.3,v=0.1)
    det  = gonio_psic.det_vectors(h=2.0,v=1.5,
                                  nu=psic.angles['nu'],
                                  delta=psic.angles['delta'])
    # get sample vectors
    sample = [[1.,1.], [.5,1.5], [-1.,1.], [-1.,-1.],[0.,.5],[1.,-1.]]
    angles = {'phi':108.0007,'chi':0.4831}
    sample = gonio_psic.sample_vectors(sample,angles=angles,gonio=psic)
    #sample = 1.5
    # compute active_area
    print active_area(psic.nm,ki=psic.ki,kr=psic.kr,beam=beam,
                      det=det,sample=sample,plot=True)

##########################################################################
if __name__ == "__main__":
    """
    test 
    """
    #test1()
    test2()
