##########################################################################
"""
Tom Trainor (tptrainor@alaska.edu)

Modifications:
--------------

"""
##########################################################################
"""
Notes

 - We should break this up into two modules:
   active_area_psic.py
   polygon_area.py

 - Make sure this handles polygons with centers off the
   origin, ie correclty handle arbitrary origin shifts

 - Interesting if we could set the detector polygon from on arbirary
   set of pixels (e.g. roi) on image detector... especially if shifted
   off center (or off center relative to a specific hkl)
"""
##########################################################################

import numpy as num

from mathutil import cosd, sind, tand
from mathutil import arccosd, arcsind, arctand
from mathutil import cartesian_mag, cartesian_angle

from gonio_psic import calc_Z, calc_D

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
      return
    elif beta < 0.0:
      print 'beta is less than 0.0'
      A_ratio = 0
      return

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
    if len(xtal) > 1:
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
        (A_beam,A_det,A_int) = _area_round_sample(beam_poly,det_poly,diameter=xtal)
    else:
        A = _area_polygon_sample(beam_poly,det_poly,sam_poly)

#########################################################################
def _area_round_sample(beam_poly,det_poly,diameter=None):
    """
    compute areas for round sample of fixed diameter
    """
    A_beam = poly_area(beam_poly)
    A_det  = poly_area(det_poly)
    # get the polygon describing the beam and detector overlap
    inner_poly = inner_polygon(beam_poly,det_poly)
    # now compute numerically the intersection area
    # of the beam on the sample within sample diameter limits
    A_int = poly_area_num(int_poly,diameter=diameter)
    return (A_beam,A_det,A_int)

#########################################################################
def _area_polygon_sample(beam_poly,det_poly,sam_poly):
    """
    compute areas for polgon sample
    """
    A_beam = poly_area(beam_poly)
    A_det  = poly_area(det_poly)
    A_sam  = poly_area(sam_poly)
    # get the polygon describing the beam and detector overlap
    inner_poly = inner_polygon(beam_poly,sam_poly)
    inner_poly = inner_polygon(inner_poly,det_poly)
    A_int = poly_area(inner_poly)
    return (A_beam,A_det,A_int)

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
def surface_intercept(k,v,diameter=None):
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
    vi = z_intercept(k,v)
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

#######################################################################
def inner_polygon(poly1,poly2):
    """
    Given two polygons find the new polygon made up from
    the inner intersections of the two
    The returned polygon is sorted.  
    """
    npts1 = len(poly1)
    npts2 = len(poly2)
    if npts1 < 3 or npts2 < 3: return None
    (poly1,angles1) = sort_points(*poly1)
    (poly2,angles2) = sort_points(*poly2)
    # loop through all possible line combinations 
    # looking for valid line intersections 
    intercepts = []
    for j in range(npts1):
        p1 = poly1[j]
        if j == npts1 - 1:
            p2 = poly1[0]
        else:
            p2 = poly1[j+1]
        for k in range(npts2):
            p3 = poly2[k]
            if k == npts2 - 1:
                p4 = poly2[0]
            else:
                p4 = poly2[k+1]
            (intercept,flag) = line_intercept(p1,p2,p3,p4)
            if flag > 0:
                intercepts.append(intercept)
    #############
    # now determine which points we can get to from
    # the origin without crossing any poly lines, 
    # ie the inner set of points
    points = []
    for p in poly1: points.append(p)
    for p in poly2: points.append(p)
    for p in intercepts: points.append(p)
    (points,angles) = sort_points(*points)
    inner_points    = []
    for p in points:
        # check against poly1
        inner = is_inner(p,poly1)
        # check against poly2
        if inner == True:
            inner = is_inner(p,poly2)
        if inner == True:
            inner_points.append(p)
    # sort the inner points
    (inner_points,angles) = sort_points(*inner_points)
    return inner_points

#######################################################################
def is_inner(point,poly):
    """
    is the point inside the polygon
    """
    npts = len(poly)
    p1 = [0.,0.]
    p2 = point
    inner = True
    k = 0
    while inner==True and k < npts:
        p3 = poly[k]
        if k == npts - 1:
            p4 = poly[0]
        else:
            p4 = poly[k+1]
        (intercept,flag) = line_intercept(p1,p2,p3,p4)
        if flag == 1 :
            inner = False
        k = k+1
    return inner

#######################################################################
def line_intercept(p1,p2,p3,p4):
    """
    Compute the point of intersection of the 2 lines defined
    by p1-p2 and p3-p4.  The p's are assumed to be 2D
    in-plane vectors (ie [x,y])

    Return ([x,y],flag).

    [x,y] = None if the lines do not intersect.  

    flag = 0 if the intercept is outside the
              extent of the 2 lines
              (or no intercept)
    flag = 1 if the intercept is in bounds,
             i.e. within the limits of the two lines
    flag = 2 if the intercept corresponds to the
             terminus of one of the lines. ie the end
             of one line is on the other line
             (inclusive of end points on end points)
    """
    # Note if vertical line m = None and b holds x-val
    (m1,b1) = line_param(p1,p2)
    (m2,b2) = line_param(p3,p4)
    if (m1 != None) and (m2 != None):
        if (m1-m2) != 0.:
            x = (b2-b1)/(m1-m2)
            y = m1*x + b1
        else:
            return (None,0)
    elif (m1 == None) and (m2 != None):
        x = b1   
        y = m2*x + b2
    elif (m1 != None) and (m2 == None):
        x = b2
        y = m1*x + b1
    else:
        return (None,0) 
    
    # min and max of points. 
    max_x1 = max(p1[0], p2[0])
    min_x1 = min(p1[0], p2[0])
    max_y1 = max(p1[1], p2[1])
    min_y1 = min(p1[1], p2[1])
    max_x2 = max(p3[0], p4[0])
    min_x2 = min(p3[0], p4[0])
    max_y2 = max(p3[1], p4[1])
    min_y2 = min(p3[1], p4[1])
    #check if the intersection is in bounds
    flag = 1
    if x > max_x1 or x < min_x1:
        flag = 0
    elif x > max_x2 or x < min_x2:
        flag = 0
    elif y > max_y1 or y < min_y1: 
        flag = 0
    elif y > max_y2 or y < min_y2: 
        flag = 0
    #check if the intersection point corresponds to an end point
    intercept = num.array([x,y])
    def _same(p1,p2,prec=0.001):
        """ are two points the same """
        #return num.all(num.equal(p1,p2))
        t1 = num.fabs(p1[0]-p2[0]) < prec
        t2 = num.fabs(p1[1]-p2[1]) < prec
        if  t1 and t2:
            return True
    if flag == 1:
        if _same(intercept,p1):
            flag = 2
        elif _same(intercept,p2):
            flag = 2
        elif _same(intercept,p3):
            flag = 2
        elif _same(intercept,p4):
            flag = 2
    return (intercept,flag)

###################################################################
def line_param(v1,v2):
    """
    Calc the params for straight line from two vectors
    (this is defined for the in-plane (2D) case)
    i.e. compute the equation of a line, y = m*x + b
    that passes through the surface point defined by v1 and v2.
    This gives back m and b.  If the line is vertical the returned
    slope (m) = None and b = x-value of the line
    """
    if (v1[0]-v2[0] != 0.):
        m = (v1[1] - v2[1])/(v1[0] - v2[0])
        b = -m*v1[0] + v1[1]
        if num.fabs(m)>1.0e6:
            m = None
            b = v1[0]
    else:  
        m =  None
        b =  v1[0]
    return (m,b)

##################################################################
def poly_area(polygon,sort=True):
    """
    Compute the area of a polygon defined
    by a list of points (xy pairs). This assumes the 
    points define a complete/enclosed polygon
    (therefore a min of 3 points).
    
    The points do not need to be in a particular order,
    this algorithm sorts them accoding to the angle
    with respect to the x-axis and then computes the
    area defined by each segment defined from each pair
    of points.

    If they are already sorted pass sort=False
    """
    npts = len(polygon)
    if npts < 3: return 0.
    if sort == True:
        (points,angles) = sort_points(*polygon)
    else:
        points = polygon
    
    # now loop through points cyclically computing
    # area of each polygon segment defined by the points
    # [0,0],[x1,y1],[x2,y2]
    A = []
    for j in range(npts):
        p1 = points[j]
        if j == npts - 1:
            p2 = points[0]
        else:
            p2 = points[j+1]
        a = segment_area(p1,p2)
        A.append(a)
    return num.sum(A)

##################################################################
def sort_points(*pts):
    """
    given a sequence of [x,y] points
    sort them according to angle (ccw w/r/t x-axis)
    return sorted list of points and angles
    """
    npts = len(pts)
    points = []
    angles = []
    #sort args by angle relative to x, c.c.w
    def _angle(v):
        # cartesian angle is always btwn 0 and 180
        angle = cartesian_angle(v,[1.,0.])
        if (v[1] < 0.):
            return 360. - angle
        else:
            return angle
    for v in pts:
        v = num.array(v[0:2])
        an = _angle(v)
        j = 0
        while j < npts -1:
            if j > len(points)-1: break
            if an < angles[j]: break
            else: j = j + 1
        points.insert(j,v)
        angles.insert(j,an)
    return (points,angles)

##################################################################
def segment_area(p1,p2):
    """
    Compute the in-plane area of the half the 
    polygon defined by the origin and two points.
    """
    # this uses cross product
    # which computes the full area of
    # the parrallogram formed by the
    # two vectors.  The polygon area
    # is half this value.  
    #p1 = num.array([p1[0],p1[1],0.])
    #p2 = num.array([p2[0],p2[1],0.])
    #a = 0.5*cartesian_mag(num.cross(p1,p2))
    # This is result of cross product
    a = (p1[0]*p2[1])**2. + (p2[0]*p1[1])**2. - (2.*p1[0]*p2[0]*p1[1]*p2[1])
    a = 0.5 * num.sqrt(a)
    return a

##################################################################
def poly_area_num(polygon,diameter=None,num_int=100,plt=False):
    """
    Numerically compute the area of a polygon
    
    If diameter != None the area is that inside the given diameter
    wrt to the center of the polygon
    
    """
    npts = len(polygon)
    if npts < 3: return 0.
    polygon = num.array(polygon)
    min_y   = min(polygon[:,1]) 
    max_y   = max(polygon[:,1]) 
    min_x   = min(polygon[:,0]) 
    max_x   = max(polygon[:,0])
    
    # compute x-values for integration
    dx = num.fabs((max_x - min_x)/float(num_int+1))
    x  = num.arange(min_x-0.5*dx, max_x+1.5*dx, dx)
    # loop through all x-vals and compute segment area
    A   = 0.0
    if plt: pline = [[],[]]
    for xx in x:
        # find polygon lines that contain x
        lines = []
        for j in range(npts):
            p1 = polygon[j]
            if j == npts - 1:
                p2 = polygon[0]
            else:
                p2 = polygon[j+1]
            if xx >= min(p1[0],p2[0]) and xx <= max(p1[0],p2[0]):
                lines.append([p1,p2])
        if len(lines) > 0:
            # now get y intercepts with vert line at xx
            p3 = [xx,min_y]
            p4 = [xx,max_y]
            y  = []
            for p1,p2 in lines:
                (inter,flag) = line_intercept(p1,p2,p3,p4)
                if flag == 0:
                    print "Error, should always get intercepts here"
                    return 0.
                else:
                    y.append(inter[1])
            numy = len(y)
            if num.mod(numy,2.) != 0 or numy == 1:
                print "Error, wrong number of intercepts!"
                return 0.
            y = num.array(y)
            y = y[num.argsort(y)]
            # now figure length of y inside the polynomial.  
            # each pair (sorted wrt y) is an inside segment.
            j = 0
            while j < numy:
                ytop = y[j+1]
                ybot = y[j]
                if diameter != None:
                    cy = (diameter/2.)**2. - xx**2.
                    if cy > 0.:
                        cy_max = num.sqrt(cy)
                        cy_min = -1.*cy_max
                        if ytop <= cy_min or ybot >= cy_max:
                            ytop = 0.0
                            ybot = 0.0
                        else:
                            if ytop > cy_max:
                                ytop = cy_max
                            if ybot < cy_min:
                                ybot = cy_min
                    else:
                        ytop = 0.0
                        ybot = 0.0
                dy = ytop - ybot
                A = A + dy*dx
                if plt:
                    pline[0].append([xx,xx])
                    pline[1].append([ybot,ytop])
                j = j+2
    # make plot of integration lines for debugging
    if plt== True:
        from matplotlib import pyplot
        for j in range(len(pline[0])):
            pyplot.plot(pline[0][j],pline[1][j],'k-')
    return A

##################################################################
def poly_y_intercepts(polygon):
    """
    find all the y intercepts of the polygon  
    """
    npts = len(polygon)
    polygon = num.array(polygon)
    min_x   = min(polygon[:,0])
    max_x   = max(polygon[:,0])
    p1      = num.array([min_x,0.])
    p2      = num.array([max_x,0.])
    intercepts = []
    for j in range(npts):
        p3 = polygon[j]
        if j == npts - 1:
            p4 = polygon[0]
        else:
            p4 = polygon[j+1]
        (intercept,flag) = line_intercept(p1,p2,p3,p4)
        if flag > 0:
            intercepts.append(intercept)
    return intercepts

##########################################################################
def plot_polygon(polygon,**kw):
    """
    plot the lines around a polygon
    """
    try:
        fmt = kw.pop('fmt')
    except:
        fmt='k'
    (points,angles) = sort_points(*polygon)
    npts = len(points)
    if npts < 3: return
    for j in range(npts):
        p1 = points[j]
        if j == npts - 1:
            p2 = points[0]
        else:
            p2 = points[j+1]
        pyplot.plot([p1[0],p2[0]],[p1[1],p2[1]],fmt,**kw)

##########################################################################
def plot_points(*pts,**kw):
    """
    plot a bunch of in-plane ([x,y]) points
    """
    try:
        fmt = kw.pop('fmt')
    except:
        fmt='k'
    npts = len(pts)
    if npts == 0: return
    xy = num.zeros((npts,2))
    for j in range(npts):
        v = pts[j]
        xy[j,0] = v[0]
        xy[j,1] = v[1]
    idx = num.argsort(xy[:,0])
    xy  = xy[idx]
    for j in range(len(xy)):
        pyplot.plot([0.,xy[j,0]],[0,xy[j,1]],fmt,**kw)

##################################################################
def plot_circle(r,**kw):
    """
    plot a circle of given radius
    """
    try:
        fmt = kw.pop('fmt')
    except:
        fmt='k'
    x = num.arange(-r,r+0.01,0.01)
    y = num.sqrt(num.fabs(r**2. - x**2.))
    pyplot.plot(x,y,fmt,**kw)
    pyplot.plot(x,-y,fmt,**kw)

##########################################################################
def trans_point(p,theta=0.,scale=1.):
    """
    simple in-plane rotation of points
    """
    M = num.array([[cosd(theta), -sind(theta)],
                   [sind(theta),  cosd(theta)]])
    pp = scale*num.dot(M,p)
    return pp

##########################################################################
##########################################################################
def test1():
    from matplotlib import pyplot
    r = 0.5
    k = [0., 1., -0.1]
    v = [.2, 0.,  1.]
    v1 = z_intercept(k,v)
    v2 = z_intercept_bounds(k,v,2*r)
    print v1,v2
    plot_circle(r)
    plot_points(v1,v2)
    pyplot.grid()

##########################################################################
def test2():
    from matplotlib import pyplot
    p1 = [1.,1.]
    p2 = [-1.,1.]
    p3 = [-1.,-1.]
    p4 = [1.,-1.]
    print poly_area([p1,p2,p3,p4])
    #
    plot_points(p1,p2,p3,p4)
    plot_polygon([p1,p2,p3,p4])
    plot_circle(2.2)
    pyplot.grid()

##########################################################################
def test3():
    from matplotlib import pyplot
    pyplot.clf()
    pyplot.grid()
    #
    poly1 = [[1.,1.], [.5,1.5], [-1.,1.], [-1.,-1.],[0.,.5],[1.,-1.]]
    poly2 = [[0.5,2.], [-0.5,2.], [-0.5,-2.],[0.5,-2.]]
    for j in range(len(poly1)):
        poly1[j] = trans_point(poly1[j],theta = -136.2,scale = .81)
    for j in range(len(poly2)):
        poly2[j] = trans_point(poly2[j],theta = 42.3,scale = 1.134)
    #
    inner = inner_polygon(poly1,poly2)
    #print inner
    #
    n = 100
    diameter=1.3
    print 'poly1 area = %6.3f, num=%6.3f' % (poly_area(poly1), poly_area_num(poly1,num_int=n))
    print 'poly2 area = %6.3f, num=%6.3f' % (poly_area(poly2), poly_area_num(poly2,num_int=n))
    print 'inner area = %6.3f, num=%6.3f' % (poly_area(inner),
                                             poly_area_num(inner,num_int=n,diameter=diameter,plt=True))
    #
    #plot_points(*poly1,**{'fmt':'ro-'})
    plot_polygon(poly1,**{'fmt':'ro-'})
    #
    #plot_points(*poly2,**{'fmt':'ko-'})
    plot_polygon(poly2,**{'fmt':'ko-'})
    #
    plot_points(*inner,**{'fmt':'go-'})
    plot_polygon(inner,fmt='g--',linewidth=4)
    #
    plot_circle(diameter/2.)

##########################################################################
if __name__ == "__main__":
    """
    test 
    """
    #test1()
    #test2()
    test3()
