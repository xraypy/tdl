##########################################################################
"""
Tom Trainor (tptrainor@alaska.edu)

Modifications:
--------------

"""
##########################################################################
"""
Notes

 - Interesting if we could set the detector polygon from on arbirary
   set of pixels (e.g. roi) on image detector... especially if shifted
   off center (or off center relative to a specific hkl)
"""
##########################################################################

import numpy as num
import types

from mathutil import cosd, sind, tand
from mathutil import arccosd, arcsind, arctand
from mathutil import cartesian_mag, cartesian_angle

from gonio_psic import calc_Z, calc_D
from polygon import poly_area, poly_area_num, inner_polygon
from polygon import plot_polygon, plot_points, plot_circle

##########################################################################
def active_area_psic_rect(gonio,slits={'wb':1.0,'hb':1.0,'wd':1.0,'hd':1.0},
                          xtal=1.,plt=False):
    """
    Calc the area of overlap of beam, sample and detector
    surface polygon projections
        A_int = area of sample illuminated within the detector projection 
        A_ratio = A_int/A_beam 
    Use to correct scattering data for area effects, including spilloff, i.e.
        Ic = I/A_ratio,   I = Idet/Io

    gonio is a gonio instance which includes the lattice parameters

    The slit settings (all in same units, eg mm):
     slits['wb'] = beam horz width 
     slits['hb'] = beam vert hieght
     slits['wd'] = detector horz width
     slits['hd'] = detector vert hieght

    If xtal is a single number we take it as the diamter of a round sample
    mounted flat.

    If xtal == None, then we assume the sample is infinite in size    
    
    Otherwise use the following structure for a general polygon shape:
    xtal = {'phi':0.,'chi':0.,'eta':0.,'mu':0.,
             'polygon':'[[],[],[],[]]'}
    polygon is a list of vectors that describe the shape of
    the xtal.  They should be given in general lab frame coordinates.
    The angle settings are the instrument angles at the time the sample
    polygon vectors were determined.

    The lab frame coordinate systems is defined such that:
        x is vertical (perpendicular, pointing to the ceiling of the hutch)
        y is directed along the incident beam path
        z make the system right handed and lies in the horizontal scattering plane
          (i.e. z is parallel to the phi axis)

    The center (0,0,0) of the lab frame is the rotation center of the instrument.  

    If the sample polygon is described at the flat phi and chi values and with
    the correct sample hieght (Z) so the sample is on the rotation center, then
    the z values of the polygon vectors will be zero.  If 2D vectors are passed
    we therefore assume these are [x,y,0]
    
    The easiest way to determine the sample coordinate vectors is to take a picture
    of the sample with a camera mounted such that is looks directly down the omega
    axis and the gonio angles set at the sample flat phi and chi values and eta = mu = 0.
    Then find the sample rotation center and measure the position of each corner (in mm)
    with up being the +x direction, and downstream being the +y direction.  

    If plt = True then makes plot

    """
    ######################
    if gonio.calc_psuedo == False:
        gonio._update_psuedo()

    alpha = gonio.pangles['alpha']
    beta = gonio.pangles['beta']

    print 'alpha = ', alpha, ', beta = ', beta
    if alpha < 0.0:
      print 'alpha is less than 0.0'
      A_ratio = 0
      return None
    elif beta < 0.0:
      print 'beta is less than 0.0'
      A_ratio = 0
      return None

    # Calc surf system transformation matrix
    M = calc_surf_transform(gonio.nm)

    #################################################################
    # transform slit vectors from lab to surface frame  
    #################################################################

    # build beam vectors, [x,y,z], in lab frame
    bh = num.array([0.,     0., 0.5*slits['wb']])
    bv = num.array([0.5*slits['hb'], 0.,  0.   ])

    # build detector vectors, [x,y,z] in lab frame
    dh = num.array([   0.,          0.,  0.5*slits['wd']])
    dv = num.array([0.5*slits['hd'],0.,        0.       ])
    # note that these assume no rotations - need to transform
    # get the rot matrix for the detector. 
    D = calc_D(nu=gonio.angles['nu'],delta=gonio.angles['delta'])
    dh = num.dot(D,dh)
    dv = num.dot(D,dv)

    # calc vectors in xs,ys,zs (i.e. surf frame)
    bh_s = num.dot(M,bh)
    bv_s = num.dot(M,bv)
    dh_s = num.dot(M,dh)
    dv_s = num.dot(M,dv)

    # create unit vectors parrallel to ki and kr 
    # and transform to surface frame.
    ki   = num.array([0.,1.,0.])
    ki_s = num.dot(M,ki)
    kr   = num.array([0.,1.,0.])
    kr   = num.dot(D,kr)
    kr_s = num.dot(M,kr)

    #################################################################
    # get the beam and detector polygons
    # calc the vector intercepts with the surface plane
    # each of the below is a vector with x_s and y_s vals 
    #################################################################
    # surface intercepts of the corners of incident beam slit 
    a = surface_intercept(ki_s,  bv_s + bh_s)
    b = surface_intercept(ki_s,  bv_s - bh_s)
    c = surface_intercept(ki_s, -bv_s - bh_s)
    d = surface_intercept(ki_s, -bv_s + bh_s)
    beam_poly = [a,b,c,d]

    # surface intercepts of the corners of the detector slit 
    e = surface_intercept(kr_s,  dv_s + dh_s)
    f = surface_intercept(kr_s,  dv_s - dh_s)
    g = surface_intercept(kr_s, -dv_s - dh_s)
    h = surface_intercept(kr_s, -dv_s + dh_s)
    det_poly = [e,f,g,h]

    #################################################################
    # get the sample position vectors
    #################################################################
    if type(xtal) == types.DictType:
        if len(xtal) < 3:
            print "Error in xtal description"
            return None
        # get sample polygon, note these are phi frame 3D vectors
        tmp = get_sample_polygon(gonio,xtal)
        sam_poly = []
        # now rotate the vectors based on the current gonio settings
        # then compute values in the surface-frame
        # Note the sample vectors SHOULD have z=0 in the surface
        # frame if they are defined correctly 
        for p_phi in tmp:
            p_m = num.dot(gonio.Z,p_phi)
            p_s = num.dot(M,p)
            if num.fabs(p_s[2]) > 0.01:
                print "Warning sample hieght problem"
            sam_poly.append(p_s[:2])
        sample_shape = True
    else:
        sample_shape = False

    #####################################################################
    # depending on whether we have a sample shape description
    # or assume a round sample, pass along appropriate
    # data to compute areas.  Note all the polygons are lists
    # of 2D surface frame vectors (ie in plane vectors)
    #####################################################################
    if shape == False:
        (A_beam,A_int) = _area_round_sample(beam_poly,det_poly,diameter=xtal)
    else:
        (A_beam,A_int) = _area_polygon_sample(beam_poly,det_poly,sam_poly)

#########################################################################
def _area_round_sample(beam_poly,det_poly,diameter=None):
    """
    compute areas for round sample of fixed diameter
    if det_poly = None, then just compute the beam and
    sample overlap ie A_int/A_beam = spill fraction
    """
    A_beam = poly_area(beam_poly)
    if det_poly != None:
        #A_det  = poly_area(det_poly)
        inner_poly = inner_polygon(beam_poly,det_poly)
        A_int = poly_area_num(inner_poly,diameter=diameter)
    else:
        A_int = poly_area_num(beam_poly,diameter=diameter)
        
    return (A_beam,A_int)

#########################################################################
def _area_polygon_sample(beam_poly,det_poly,sam_poly):
    """
    compute areas for polgon sample
    if det_poly = None, then just compute the beam and
    sample overlap ie A_int/A_beam = spill fraction
    """
    A_beam = poly_area(beam_poly)
    #A_sam  = poly_area(sam_poly)
    inner_poly = inner_polygon(beam_poly,sam_poly)
    if det_poly != None:
        #A_det  = poly_area(det_poly)
        inner_poly = inner_polygon(inner_poly,det_poly)
    A_int = poly_area(inner_poly)
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

    Note that nm is the surface normal (unit vector), defined
    in the laboratory m-frame (i.e. pointing in an arbitrary
    direction for a particular gonio setting) - see gonio_psic.py
    for more details.

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
    v_z = nm

    v   = num.array([0.,-1.,0.])
    v_y = v - (num.dot(v,v_z) * v_z) / (cartesian_mag(v_z)**2.) 
    v_y = v_y / cartesian_mag(v_y)
    
    v_x = num.cross(v_y,v_z)
    v_x = v_x / cartesian_mag(v_x)

    F = num.array([v_x,v_y,v_z])
    M = num.linalg.inv(F.transpose())

    return M

##########################################################################
def get_sample_polygon(xtal):
    """
    xtal = {'phi':0.,'chi':0.,'eta':0.,'mu':0.,
             'polygon':'[[],[],[],[]]'}
    polygon is a list of vectors that describe the shape of
    the xtal.  They should be given in general lab frame coordinates.
    The angle settings are the instrument angles at the time the sample
    polygon vectors were determined.
    
    The lab frame coordinate systems is defined such that:
        x is vertical (perpendicular, pointing to the ceiling of the hutch)
        y is directed along the incident beam path
        z make the system right handed and lies in the horizontal scattering plane
          (i.e. z is parallel to the phi axis)

    The center (0,0,0) of the lab frame is the rotation center of the instrument.

    If the sample polygon is described at the flat phi and chi values and with
    the correct sample hieght (Z) so the sample is on the rotation center, then
    the z values of the polygon vectors will be zero.  If 2D vectors are passed
    we therefore assume these are [x,y,0]

    Note that the sample_poly that is returned is a list of 3D vectors
    defined in the lab phi frame.  These will need to be rotated to the
    M-frame for an arbitrary gonio setting
    """
    # get sample rotation matrix
    Z = calc_Z(phi=xtal['phi'],chi=xtal['chi'],eta=xtal['eta'],mu=xtal['mu'])
    Zinv = num.linalg.inv(Z)

    polygon_phi = []
    polygon = xtal['polygon']
    if len(polygon) < 3:
        print "Sample polygon must be 3 or more points"
        return None
    
    # If p's have anly two components then we assume they are xy pairs
    # therefore we can add a third value of zero for z
    # Then since p's are defined in lab frame at given set of angles
    # we need to unrotate to get back to phi frame vectors
    for p in polygon:
        if len(p) == 2: p = p.resize((3,))
        p_phi = num.dot(Zinv,p)
        polygon_phi.append(p_phi)
        
    return polygon_phi

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
    from matplotlib import pyplot
    r = 0.5
    k = [0., 1., -0.1]
    v = [.2, 0.,  1.]
    v1 = surface_intercept(k,v)
    v2 = surface_intercept_bounds(k,v,2*r)
    print v1,v2
    plot_circle(r)
    plot_points(v1,v2)
    pyplot.grid()

##########################################################################
if __name__ == "__main__":
    """
    test 
    """
    test1()

