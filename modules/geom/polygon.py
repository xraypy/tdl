"""
Generalized Polygon Computations

Authors/Modifications:
----------------------
* Tom Trainor (tptrainor@alaska.edu)

Todo:
-----
* Make sure this handles polygons with centers off the
  origin, ie correclty handle arbitrary origin shifts

* Add some more polygon calcs
  - compute center
  - determine type (simple, complex etc)
"""
##########################################################################

import numpy as num
from matplotlib import pyplot

from tdl.modules.utils.mathutil import cosd, sind, tand
from tdl.modules.utils.mathutil import arccosd, arcsind, arctand
from tdl.modules.utils.mathutil import cartesian_mag, cartesian_angle

#########################################################################
def inner_polygon(poly1,poly2):
    """
    Given two polygons find the new polygon made up from
    the inner intersections of the two

    Parameters:
    -----------
    * poly1 and poly2 are polygons defined by 3 or more
      2D vectors (ie points in an xy plane)

    Outputs:
    --------
    * The returned polygon is sorted.

    Example:
    --------
    >>poly1 = [[1.,1.], [.5,1.5], [-1.,1.], [-1.,-1.],[0.,.5],[1.,-1.]]
    >>poly2 = [[0.5,2.], [-0.5,2.], [-0.5,-2.],[0.5,-2.]]
    >>inner = polygon.inner_polygon(poly1,poly2)
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
    Is the point inside the polygon

    Parameters:
    -----------
    * point is a [x,y] pair
    * poly is a list of 3 or more [x,y] pairs defining a polygon
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
    by p1-p2 and p3-p4.

    Parameters:
    -----------
    * The p's are 2D in-plane vectors (ie [x,y])

    Output: ([x,y],flag)
    -------
    * [x,y] is the point of intersection or
      None if the lines do not intersect.  
    * flag = 0 if the intercept is outside the
      extent of the 2 lines (or no intercept)
    * flag = 1 if the intercept is in bounds,
      i.e. within the limits of the two lines
    * flag = 2 if the intercept corresponds to the
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
    def _same(p1,p2,prec=0.00001):
        """ are two points the same """
        #return num.all(num.equal(p1,p2))
        t1 = num.fabs(p1[0]-p2[0]) < prec
        t2 = num.fabs(p1[1]-p2[1]) < prec
        if  t1 and t2:
            #print "same", p1,p2
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
    Calc the params for straight line defined from two vectors
    (this is defined for the in-plane (2D) case)

    Parameters:
    -----------
    * v1 and v2 are 2D vectors ([x,y])

    Outputs: (m,b)
    --------
    * m is the line slope
    * b is the line y-intercept

    Notes:
    ------
    This computes the equation of a line, y = m*x + b,
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
    Compute the area of a polygon

    Parameters:
    -----------
    * polygon is a list of points (xy pairs).
      This assumes the points define a complete/
      enclosed polygon (therefore a min of 3 points).

    * sort is a flag to indictate if the points
      should be sorted by thier angles relative to
      the x-axis

    Notes:
    ------
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
    Sort points according to angle (ccw w/r/t x-axis)

    Parameters:
    -----------
    * a sequence of [x,y] points

    Outputs: (points, angles)
    --------
    * return sorted list of points and angles
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
    Compute the in-plane area of the polygon defined
    by the origin and two points (p1 and p2).
    """
    # this uses cross product
    # which computes the full area of
    # the parrallogram formed by the
    # two vectors.  The polygon area
    # is half this value.  
    #p1 = num.array([p1[0],p1[1],0.])
    #p2 = num.array([p2[0],p2[1],0.])
    #a = 0.5*cartesian_mag(num.cross(p1,p2))
    # This is the result of the cross product operation:
    a = (p1[0]*p2[1])**2. + (p2[0]*p1[1])**2. - (2.*p1[0]*p2[0]*p1[1]*p2[1])
    a = 0.5 * num.sqrt(a)
    return a

##################################################################
def poly_area_num(polygon,diameter=None,num_int=100,plot=False):
    """
    Numerically compute the area of a polygon

    Parameters:
    -----------
    * polygon is a list of [x,y] points defining the shape
    * if diameter != None the area is that inside the given diameter
      wrt to the center of the polygon
    * num_int is the number of divisions to use in the integration
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
    if plot: pline = [[],[]]
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
                if plot:
                    pline[0].append([xx,xx])
                    pline[1].append([ybot,ytop])
                j = j+2
    # make plot of integration lines for debugging
    if plot== True:
        for j in range(len(pline[0])):
            pyplot.plot(pline[0][j],pline[1][j],'k-')
    return A

##################################################################
def poly_y_intercepts(polygon):
    """
    find all the y-axis intercepts of the polygon  
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
def trans_point(p,theta=0.,scale=1.):
    """
    simple in-plane rotation of points
    """
    M = num.array([[cosd(theta), -sind(theta)],
                   [sind(theta),  cosd(theta)]])
    pp = scale*num.dot(M,p)
    return pp

##########################################################################
def plot_polygon(polygon,**kw):
    """
    plot the lines around a polygon
    """
    try:
        fmt = kw.pop('fmt')
    except:
        fmt='k'
    try:
        label = kw.pop('label')
    except:
        label = None
    (points,angles) = sort_points(*polygon)
    npts = len(points)
    if npts < 3: return
    for j in range(npts):
        p1 = points[j]
        if j == npts - 1:
            p2 = points[0]
        else:
            p2 = points[j+1]
        if j < npts - 1:
            pyplot.plot([p1[0],p2[0]],[p1[1],p2[1]],fmt,**kw)
        else:
            pyplot.plot([p1[0],p2[0]],[p1[1],p2[1]],fmt,label=label,**kw)

##########################################################################
def plot_points(points,**kw):
    """
    plot a bunch of in-plane ([x,y]) points

    Parameters:
    -----------
    * points is a list or tupe of [x,y] pairs
    """
    try:
        fmt = kw.pop('fmt')
    except:
        fmt='k'
    try:
        label = kw.pop('label')
    except:
        label = None
    npts = len(points)
    if npts == 0: return
    xy = num.zeros((npts,2))
    for j in range(npts):
        v = points[j]
        xy[j,0] = v[0]
        xy[j,1] = v[1]
    idx = num.argsort(xy[:,0])
    xy  = xy[idx]
    for j in range(len(xy)):
        if j < npts - 1:
            pyplot.plot([0.,xy[j,0]],[0,xy[j,1]],fmt,**kw)
        else:
            pyplot.plot([0.,xy[j,0]],[0,xy[j,1]],fmt,label=label,**kw)

##################################################################
def plot_circle(r,**kw):
    """
    plot a circle of given radius
    """
    try:
        fmt = kw.pop('fmt')
    except:
        fmt='k'
    try:
        label = kw.pop('label')
    except:
        label = None
    x = num.arange(-r,r+0.01,0.01)
    y = num.sqrt(num.fabs(r**2. - x**2.))
    pyplot.plot(x,y,fmt,**kw)
    pyplot.plot(x,-y,fmt,label=label,**kw)

##########################################################################
##########################################################################
def test():
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
    print 'poly1 area = %6.3f, num=%6.3f' % (poly_area(poly1),
                                             poly_area_num(poly1,num_int=n))
    print 'poly2 area = %6.3f, num=%6.3f' % (poly_area(poly2),
                                             poly_area_num(poly2,num_int=n))
    print 'inner area = %6.3f, num=%6.3f' % (poly_area(inner),
                                             poly_area_num(inner,num_int=n,
                                             diameter=diameter,plot=True))
    #
    plot_polygon(poly1,fmt='ro-')
    plot_polygon(poly2,fmt='ko-')
    plot_points(inner,fmt='go-')
    plot_polygon(inner,fmt='g--',linewidth=4)
    plot_circle(diameter/2.)

##########################################################################
if __name__ == "__main__":
    """
    test 
    """
    test()
