##########################################################################
"""
Tom Trainor (tptrainor@alaska.edu)

Modifications:
--------------

"""
##########################################################################
"""
Notes

"""
##########################################################################

import numpy as num

from mathutil import cosd, sind, tand
from mathutil import arccosd, arcsind, arctand
from mathutil import cartesian_mag, cartesian_angle

##########################################################################

def calc_surf_transform(nm):
    """
    This routine calculates the matrix, M, which transforms the
    indicies of vectors defined in the lab frame basis to the
    a surface frame defined such that:
        zs is along the surface normal (parrallel to nm)
        ys is the projection of the -y axis onto the surface
        xs is defined to make it a right handed orthonormal set
           (and is therefore in the surface plane)

    Note that nm is the surface normal (unit vector) pointing in
    an arbitrary definition defined by a particular gonio setting.  
    The vector nm is defined in the psic laboratory basis.

    The transform is defined as:
        |e_xs|      |e_x|              |vx|     
        |e_ys| = F* |e_y|    and   F = |vy|   
        |e_zs|      |e_z|              |vz|    
    Where e's are the (cartesian) basis vectors.  Given a
    vector [x,y,z] defined in the lab fram basis [e_x,e_y,e_z],
    we can compute the indicies of this vector in the surface
    basis from:
        |xs|      |x|
        |ys| = M* |y|
        |zs|      |z|
    where
        M = transpose(inv(F)) = inv(transpose(F))
        
    Note see lattice.py for more notes on general transoform, and see
    gonio_psic.py for note on calc of the surface normal vector.  
    """
    v_z = gonio.nm

    v   = num.array([0.,-1.,0.])
    v_y = v - (num.dot(v,v_z) * v_z) / (cartesian_mag(v_z)**2.) 
    v_y = v_y / cartesian_mag(v_y)
    
    v_x = num.cross(v_y,v_z)
    v_x = v_x / cartesian_mag(v_x)

    F = num.array([v_x,v_y,v_z])
    M = num.linalg.inv(F.transpose())

    return M

##########################################################################
def active_area_psic_rect(gonio,p_phi,p_chi,wb,hb,wd,hd,xtal_shape,plt);
#def active_area_psic_rect(angles,Uarray,n_vec,n_flag,p_phi,p_chi,wb,hb,wd,hd,xtal_shape,plt);
    """
    Calc the area of overlap of beam, sample and det parrallelograms
        A_ratio = intersection_area/beam_area 
    Use to correct scattering data for area effects, i.e.
        Ic = I/A_ratio,   I = Idet/Io

    gonio is a gonio instance which includes the lattice parameters

    the slit settings (all in same units, eg mm):
    wb = beam horz width 
    hb = beam vert hieght
    wd = detector horz width
    hd = detector vert hieght

    xtal_shape = matrix giving co-ordinates of the four corners of the sample 
                 at all angle =0 and phi=f_phi and chi= f_chi in lab-frame cartesian
                 co-ordinate ( as beam direction is +ve y, up- is +ve x and outboard is z)
      
    if plt = 1 then makes plot


    see: m_files\sanjit\area_psic_rect_num.m
         m_files\current\saag\xrd_ctr\data_red\int_ctr_data_3\parse_int_mac_rect.m
    """
    ######################
    if gonio.calc_psuedo == False:
        gonio._update_psuedo()

    alpha = gonio.pangles['alpha']
    beta = gonio.pangles['beta']

    print 'alpha = ', alpha, ', beta = ', beta
    if alpha < -0.0:
      print 'alpha is less than 0.0'
      A_ratio = 0
      return
    elif beta < -0.0:
      print 'beta is less than -0.0'
      A_ratio = 0
      return

    # Calc surf sys. transformation matrix
    M = calc_surf_transform(gonio.nm)

    ###############################################################################
    ####   NEW ALOGORTHM INCLUDED FOR A FOUR SIDED SAMPLE SHAPE #########
    #################################################################################
    # Sample positions determined from the picture taken at 
    #del=0;eta=0;chi=f_chi;phi=f_phi;nu= 0;mu=0
    ### xtal_shape =[ x1 x2 x3 x4     sample four corner co-ordinates in cartesian co-ordinate (lab-frame)
    #                 y1 y2 y3 y4     determined from the camera picture and size
    #                 z1 z2 z3 z4]    assuming center to be (0 0 0)
    #
    ######### If sample co-ordinates are defined at f_phi & f_chi
    # rotate the vector from f_phi & f_chi position to zero lab-frame position
    # with a new Z matrix with following angular values
    #del=0;eta=0;chi=f_chi;phi=f_phi;nu= 0;mu=0
    angles_0 = [0,0,p_chi,p_phi,0,0];
    k = 1;
    [Z_0,Qm_0,ki_0,kr_0,tth_0] = calc_Z_and_Q(angles_0,k);
    #Z_0

    
    xtal_lab_unrot = xtal_shape';
    Svc_1 = xtal_lab_unrot(1:3,1);
    Svc_2 = xtal_lab_unrot(4:6,1);
    Svc_3 = xtal_lab_unrot(7:9,1);
    Svc_4 = xtal_lab_unrot(10:12,1);
    # get the zero lab frame positions for the sample vectors
    Svc_1 = Z_0*Svc_1;
    Svc_2 = Z_0*Svc_2;
    Svc_3 = Z_0*Svc_3;
    Svc_4 = Z_0*Svc_4;
    #xtal__lab0_coord = [ Svc_1 Svc_2 Svc_3 Svc_4]

    ### after rotation of the crystal: rotated base
    Svcr_1 = Z*Svc_1;
    Svcr_2 = Z*Svc_2;
    Svcr_3 = Z*Svc_3;
    Svcr_4 = Z*Svc_4;

    #### vectors transformed in surf-frame
    Svs_1 = M*Svcr_1;
    Svs_2 = M*Svcr_2;
    Svs_3 = M*Svcr_3;
    Svs_4 = M*Svcr_4;

    #xtal_in_surf_coord = [ Svs_1 Svs_2 Svs_3 Svs_4 ]

    #################################################################
    # transform vectors from lab to surface frame                   #
    #################################################################

    # build beam vectors in x,y,z (i.e. lab frame)
    bh = [  0    ;   0;   0.5*wb];
    bv = [0.5*hb ;   0;     0   ];

    # build detector vectors in x,y,z (i.e. lab frame)
    # note that these assume no rotations, need to transform as below
    dh = [       0   ;   0 ;   0.5*wd];
    dv = [    0.5*hd ;   0 ;     0   ];

    # calc the rot matrix for the detector. 
    del = angles(1);
    nu = angles(5);
    DM_1 = [  cosd(del),  sind(del),  0  ; 
             -sind(del),  cosd(del),  0  ;
             0     ,     0     ,  1 ];
          
    DM_2 = [    1 ,     0   ,      0     ;
                0,  cosd(nu), -sind(nu)  ; 
                0,  sind(nu),  cosd(nu) ];
          
    DM = (DM_2)*(DM_1);

    dh = DM*dh;
    dv = DM*dv;

    # calc vectors in xs, ys,zs (i.e. surf frame)
    bh_s = M*bh;
    bv_s = M*bv;

    dh_s = M*dh;
    dv_s = M*dv;
    #disp('******* CHECKED UP TO HERE for syntax,dimensions and  results *******')

    # note ki and kr calc above and are already unit vectors
    # now transform to surface frame.
    ki_s = M*ki;
    kr_s = M*kr;

    ###################################################################
    # calc the intercepts of the corners of the beam,detector and sample        #
    # each of the below is a row vector with x and y vals             #
    ###################################################################

    [a] = z_intercept(ki_s, [ bv_s + bh_s]);
    [b] = z_intercept(ki_s, [ bv_s - bh_s]);
    [c] = z_intercept(ki_s, [-bv_s - bh_s]);
    [d] = z_intercept(ki_s, [-bv_s + bh_s]);

    # calc the intercepts of the corners of the detector 

    [e] = z_intercept(kr_s, [ dv_s + dh_s]);
    [f] = z_intercept(kr_s, [ dv_s - dh_s]);
    [g] = z_intercept(kr_s, [-dv_s - dh_s]);
    [h] = z_intercept(kr_s, [-dv_s + dh_s]);

    # calc the intercepts of the corners of the sample ###
    # make these two dimension for cross product below
    i= [Svs_1(1,1), Svs_1(2,1)];
    j= [Svs_2(1,1), Svs_2(2,1)];
    k= [Svs_3(1,1), Svs_3(2,1)];
    l= [Svs_4(1,1), Svs_4(2,1)];

    #vortices = [a;b;c;d;e;f;g;h;i;j;k;l];

    ############## NUM CALCULATION OF INTERSECTION AREA ##########################
    #Use the minimum intersection area formed out of beam,sample and detector ###
    #the intersection area is defined by the ab,bc,cd,da,ef,fg,gh,he,ij,jk,kl & #
    #lj lines #
    #############################################################################

    # params for lines                                                 #
    # the lines array has a row for each intercept with                #
    # col1 = m, col2 = b.  Note if a slope is inf or 0. large          #
    # m is returned as inf and the x-val for the line is in b          #

    ####################################################################
    #### if any line parallel to x-axis then (slope m = 0)
    #### check for slope <= 0.1, one can change it to lower/higher value
    #### rotate all the vertices by a rotation angle

    lines =line_foot_print(a,b,c,d,e,f,g,h,i,j,k,l);
    l_idx = length(find(abs(lines(:,1))<=0.1));
    m_rot_angle = 1;
    while l_idx ~=0, 
       m_rot_mat= [cosd(m_rot_angle),-sind(m_rot_angle);...
                   sind(m_rot_angle),cosd(m_rot_angle)];
        
       a = (m_rot_mat*a(:))';
       b = (m_rot_mat*b(:))';
       c = (m_rot_mat*c(:))';
       d = (m_rot_mat*d(:))';
       e = (m_rot_mat*e(:))';
       f = (m_rot_mat*f(:))';
       g = (m_rot_mat*g(:))';
       h = (m_rot_mat*h(:))';
       i = (m_rot_mat*i(:))';
       j = (m_rot_mat*j(:))';
       k = (m_rot_mat*k(:))';
       l = (m_rot_mat*l(:))';
       lines= line_foot_print(a,b,c,d,e,f,g,h,i,j,k,l);
       l_idx = length(find(abs(lines(:,1))<=0.1));
    end

    #vortices = [a;b;c;d;e;f;g;h;i;j;k;l];
    lines =line_foot_print(a,b,c,d,e,f,g,h,i,j,k,l);

    ####################################################################
    # calc the area of the beam and det print                          #
    # note these use long vals ie. tot f_print including beyond xtal   #
    ####################################################################
    # here just calc areas
    #disp('******* beam,detector and sample footprint 2D area *******')

    beam_area_surf =  norm(cross([a,0],[b,0])) ...
                     + norm(cross([b,0],[c,0]));
    det_area_surf  = norm(cross([e,0],[f,0])) ...
                     + norm(cross([f,0],[g,0]));
    sample_area_surf = norm(cross([i,0],[j,0])) ...
                     + norm(cross([k,0],[l,0]));

    ###############################################################
    # find the y-intercepts( x-co-ordinate) for  all the lines
    # until the x-height falls below tollerance
    # then integrate the area under each del_y, area=del_y*x_height
    # do the intergartion from origin to +ve y values first
    # then origin to negavtive y values and then add the two to get the final
    # area of the intersection polygon
    #####################################################

    ###############Integration for the positive side from the center(0,0)##

    start_ints_area = min([beam_area_surf,det_area_surf,sample_area_surf]);
    num_steps_p = 100; # number of steps one side
    del_y = (sqrt(start_ints_area)/2)/(num_steps_p+1); # step size
    y_ints=0.0;
    line_yints=line_y_intercepts(y_ints,lines);
    l_line_yints = length(line_yints);
    sortx_line_yints= sort(line_yints(:,1));
    x_height= abs(sortx_line_yints(7)- sortx_line_yints(6));
    sum_area_p=y_ints*x_height;
    count_p=0;
    count_p_s=0;

    while x_height >= 0.001, # 0.001 is the tollerance of the x-height
       y_ints=y_ints+del_y;
       line_yints=line_y_intercepts(y_ints,lines);
       l_line_yints = length(line_yints);
       sortx_line_yints= sort(line_yints(:,1));
       x_height_1 = x_height;
       x_height= abs(sortx_line_yints(7)- sortx_line_yints(6));
       sum_area_p = sum_area_p + (del_y*x_height);
       # rescaling the slices
       if (x_height_1 - x_height)>= 0.1*x_height_1;
          del_y=del_y*0.2;
          count_p_s=count_p_s+1;
       end
      count_p=count_p+1;
    end
    #count_p
    #x_height_f=x_height
    #y_ints_f=y_ints
    #sum_area_p_f=sum_area_p
    #sum_area_p
    #######################################################################

    ###############Integration for the negative side from the center(0,0)##

    del_y = (sqrt(start_ints_area)/2)/(num_steps_p+1); # step size
    y_ints=0.0;
    line_yints=line_y_intercepts(y_ints,lines);
    l_line_yints = length(line_yints);
    sortx_line_yints= sort(line_yints(:,1));
    x_height= abs(sortx_line_yints(7)- sortx_line_yints(6));
    sum_area_n=y_ints*x_height;
    count_n=0;
    count_n_s=0;

    while x_height >= 0.001,# 0.001 is the tollerance of the x-height
       y_ints=y_ints-del_y;
       line_yints=line_y_intercepts(y_ints,lines);
       l_line_yints = length(line_yints);
       sortx_line_yints= sort(line_yints(:,1));
       x_height_1 = x_height;
       x_height= abs(sortx_line_yints(7)- sortx_line_yints(6));
       sum_area_n = sum_area_n + (abs(del_y)*x_height);
       # rescaling the slices
       if (x_height_1 - x_height)>= 0.1*x_height_1;
          del_y=del_y*0.2;
          count_n_s=count_n_s+1;
       end
       count_n = count_n+1;
    end 
    #count_n
    #sum_area_n
    bsd_ints_area = sum_area_p + sum_area_n;
    A_ratio = bsd_ints_area/beam_area_surf;

    # include the areas of the sample, beam and detector slits
    t= ['Area correction = ' num2str(A_ratio)];
    disp(t)

    if plt ==1
     plot_lines(a,b,c,d,e,f,g,h,i,j,k,l);
    end
      
    ###############get the boundary of integration ##################
    #y_corner = [a(:,2),b(:,2),c(:,2),d(:,2),e(:,2),f(:,2),...
    #               g(:,2),h(:,2),i(:,2),j(:,2),k(:,2),l(:,2)];
    #[y_pos_min,y_pos_max,y_neg_min,y_neg_min] = sort_matrix(y_corner);

    ###############################################################################

     
#################################################################
################## functions ####################################
#################################################################

#find min max of a matrix not including zero
function [pos_min,pos_max,neg_min,neg_max] = sort_matrix(matrix)
    matrix = sort(matrix);
    pos_idx = find(matrix(:)>0);
    l_posidx = length(find(matrix(:)>0));           
    pos_min = min(matrix(:,pos_idx(1,1):pos_idx(l_posidx,1)));
    pos_max = max(matrix(:,pos_idx(1,1):pos_idx(l_posidx,1)));

    neg_idx = find(matrix(:)<0);
    l_negidx = length(find(matrix(:)<0));          
    neg_min = min(matrix(:,neg_idx(1,1):neg_idx(l_negidx,1)));
    neg_max = max(matrix(:,neg_idx(1,1):neg_idx(l_negidx,1)));

###########################################################
function line_yints = line_y_intercepts(y_ints,lines)
    for i=1:1:12
    line1=[0,y_ints];
    line2=lines(i,:);
    [line_yints(i,1),line_yints(i,2)] = y_line_intercept( line1,line2);
    end

    ##################################################################
    # find the x-value  between lines y = b1(parallel to x) and y = m2*x +b2 
    function [x,y] = y_line_intercept( line1,line2)

    # m1=0 and b1 are params for line parallel to x
    # m2 and b2 are params for line2
    # therefore find the intersections of the line and see if 
    # and get the x-value

    m1 = line1(1);
    b1 = line1(2);
    m2 = line2(1);
    b2 = line2(2);


    if (m1~=inf) & (m2~=inf);
    if (m1-m2)~=0
     x= (b2-b1)/(m1-m2);
     y= m1*x + b1;
    else
     # x= (b2-b1)/eps;
     # y= m1*x + b1;
     x= inf;
     y = inf;
    end;
    elseif (m1==inf) & (m2~=inf)
    x = b1;   # recall if vertical line b holds x-val
    y = m2*x + b2;
    elseif (m1~=inf) & (m2==inf)
    x = b2;   # recall if vertical line b holds x-val
    y = m1*x + b1;
    else;
    x = inf;
    y = inf;
    end;

##################################################################
function [int] = z_intercept(k,v)

    # find x,y for z=0
    ks = [k(1),k(2)];
    x_int = v(1) - (k(1)/k(3))*v(3);
    y_int = v(2) - (k(2)/k(3))*v(3);

    int = [x_int, y_int];

##################################################################
# black sample
# red is beam
# green is det
function  plot_lines(a,b,c,d,e,f,g,h,i,j,k,l)                           
    disp('black is sample')
    disp('red is beam')
    disp('green is det')

    clf, hold on, gup

    plot([a(1),b(1)],[a(2), b(2)],'r')
    plot([b(1),c(1)],[b(2), c(2)],'r')
    plot([c(1),d(1)],[c(2), d(2)],'r')
    plot([d(1),a(1)],[d(2), a(2)],'r')

    plot([e(1),f(1)],[e(2), f(2)],'g')
    plot([f(1),g(1)],[f(2), g(2)],'g')
    plot([g(1),h(1)],[g(2), h(2)],'g')
    plot([h(1),e(1)],[h(2), e(2)],'g')

    plot([i(1),j(1)],[i(2), j(2)],'k')
    plot([j(1),k(1)],[j(2), k(2)],'k')
    plot([k(1),l(1)],[k(2), l(2)],'k')
    plot([l(1),i(1)],[l(2), i(2)],'k')

###################################################################
function [m,b] = line_param(v1,v2);

    # calc params for straight line from two vectors
    # i.e. y = m*x + b , this gives back m and b

    if (v1(1)-v2(1) ~= 0)
    m = (v1(2) - v2(2))/(v1(1) - v2(1));
    b = -m*v1(1) + v1(2);
    if abs(m)>1e6;     # arbitrary large number!!!!
     m = inf;         # put these here to deal with numerical problems
     b = v1(1);       # when you have near vert. lines (i.e. make em vert)
    # disp('inf1')
    end; 
    else;  # set m to inf and b to the x-val for vertical line!!!
    m =  inf;
    b =  v1(1);
    # disp('inf2')
    end;

#######################################################################
function line_ints = get_line_intercepts(lines, a,b,c,d,e,f,g,h)

    I(1,:) = line_intercept(lines(1,:), a, b, lines(5,:), e, f);
    I(2,:) = line_intercept(lines(1,:), a, b, lines(6,:), f, g);
    I(3,:) = line_intercept(lines(1,:), a, b, lines(7,:), g, h);
    I(4,:) = line_intercept(lines(1,:), a, b, lines(8,:), h, e);
    I(5,:) = line_intercept(lines(2,:), b, c, lines(5,:), e, f);
    I(6,:) = line_intercept(lines(2,:), b, c, lines(6,:), f, g);
    I(7,:) = line_intercept(lines(2,:), b, c, lines(7,:), g, h);
    I(8,:) = line_intercept(lines(2,:), b, c, lines(8,:), h, e);
    I(9,:) = line_intercept(lines(3,:), c, d, lines(5,:), e, f);
    I(10,:) = line_intercept(lines(3,:), c, d, lines(6,:), f, g);
    I(11,:) = line_intercept(lines(3,:), c, d, lines(7,:), g, h);
    I(12,:) = line_intercept(lines(3,:), c, d, lines(8,:), h, e);
    I(13,:) = line_intercept(lines(4,:), d, a, lines(5,:), e, f);
    I(14,:) = line_intercept(lines(4,:), d, a, lines(6,:), f, g);
    I(15,:) = line_intercept(lines(4,:), d, a, lines(7,:), g, h);
    I(16,:) = line_intercept(lines(4,:), d, a, lines(8,:), h, e);

    line_ints = I;


#######################################################################
function int = line_intercept( line1, v1_1, v1_2, line2, v2_1, v2_2);

    # m1 and b1 are params for line joining points v1_1 and v1_2
    # m2 and b2 are params for line joining points v2_1 and v2_2
    # therefore find the intersections of the line and see if 
    # the intersection is within the bounds of the lines

    # int is returned with a int(3) = 0 if the intersection
    # is in bounds, i.e. within the limits of the two lines

    m1 = line1(1);
    b1 = line1(2);
    m2 = line2(1);
    b2 = line2(2);


    if (m1~=inf) & (m2~=inf);
    if (m1-m2)~=0
     x= (b2-b1)/(m1-m2);
     y= m1*x + b1;
    else
     # x= (b2-b1)/eps;
     # y= m1*x + b1;
     x= inf;
     y = inf;
    end;
    elseif (m1==inf) & (m2~=inf)
    x = b1;   # recall if vertical line b holds x-val
    y = m2*x + b2;
    elseif (m1~=inf) & (m2==inf)
    x = b2;   # recall if vertical line b holds x-val
    y = m1*x + b1;
    else;
    x = inf;
    y = inf;
    end;

    # limits of line 2
    hx2 = max(v2_1(1), v2_2(1));
    lx2 = min(v2_1(1), v2_2(1));
    hy2 = max(v2_1(2), v2_2(2));
    ly2 = min(v2_1(2), v2_2(2));

    # limits of line 1
    hx1 = max(v1_1(1), v1_2(1));
    lx1 = min(v1_1(1), v1_2(1));
    hy1 = max(v1_1(2), v1_2(2));
    ly1 = min(v1_1(2), v1_2(2));

    # check if the intersection is outside the range of the lines 
    #if (x > hx2) | (x < lx2) | (x > hx1) | (x < lx1)
    #  int = [x,y,1];
    #elseif  (y > hy2) | (y < ly2) | (y > hy1) | (y < ly1) 
    #  int = [x,y,1];
    #else;
    #  int = [x,y,0];
    #end

    p=0.001;

    if cmp(x,hx2,p)==1 | cmp(x,lx2,p)==-1 | cmp(x,hx1,p)==1 | cmp(x,lx1,p)==-1
    int = [x,y,1];
    elseif cmp(y,hy2,p)==1 | cmp(y,ly2,p)==-1 | cmp(y,hy1,p)==1 | cmp(y,ly1,p)==-1 
    int = [x,y,1];
    else;
    int = [x,y,0];
    end



############################################################################
function int_1 = line_intercept_1( line1, v1_1, v1_2, line2, v2_1, v2_2);

    # m1 and b1 are params for line joining points v1_1 and v1_2
    # m2 and b2 are params for line joining points v2_1 and v2_2
    # therefore find the intersections of the line and see if 
    # the intersection is within the bounds of the lines

    # int is returned with a int(3) = 0 if the intersection
    # is in bounds, i.e. within the limits of the two lines

    m1 = line1(1);
    b1 = line1(2);
    m2 = line2(1);
    b2 = line2(2);


    if (m1~=inf) & (m2~=inf);
    if (m1-m2)~=0
     x= (b2-b1)/(m1-m2);
     y= m1*x + b1;
    else
     # x= (b2-b1)/eps;
     # y= m1*x + b1;
     x= inf;
     y = inf;
    end;
    elseif (m1==inf) & (m2~=inf)
    x = b1;   # recall if vertical line b holds x-val
    y = m2*x + b2;
    elseif (m1~=inf) & (m2==inf)
    x = b2;   # recall if vertical line b holds x-val
    y = m1*x + b1;
    else;
    x = inf;
    y = inf;
    end;

    # limits of line 2
    hx2 = max(v2_1(1), v2_2(1));
    lx2 = min(v2_1(1), v2_2(1));
    hy2 = max(v2_1(2), v2_2(2));
    ly2 = min(v2_1(2), v2_2(2));

    # limits of line 1
    hx1 = max(v1_1(1), v1_2(1));
    lx1 = min(v1_1(1), v1_2(1));
    hy1 = max(v1_1(2), v1_2(2));
    ly1 = min(v1_1(2), v1_2(2));

    # check if the intersection is outside the range of the lines 
    #if (x > hx2) | (x < lx2) | (x > hx1) | (x < lx1)
    #  int = [x,y,1];
    #elseif  (y > hy2) | (y < ly2) | (y > hy1) | (y < ly1) 
    #  int = [x,y,1];
    #else;
    #  int = [x,y,0];
    #end

    p=0.001;

    if cmp(x,hx2,p)==1 | cmp(x,lx2,p)==-1 | cmp(x,hx1,p)==1 | cmp(x,lx1,p)==-1
     int_1 = [x,y,1];
    elseif cmp(y,hy2,p)==1 | cmp(y,ly2,p)==-1 | cmp(y,hy1,p)==1 | cmp(y,ly1,p)==-1 
    int_1 = [x,y,1];
    else;
    int_1 = [x,y,0];
    end



############################################################################
function c = cmp(a,b,prec)
    # compare two numbers
    # c = 0 if a = b within prec
    # c = 1 if a > b
    # c = -1 if a < b
    # maybe problems with infs

    if abs(a-b) <= prec
    c = 0;
    elseif a > b;
    c = 1;
    elseif a < b;
    c = -1;
    end;

#############################################################################
function y_ints = get_y_ints(lines,a,b,c,d,e,f,g,h)

    yz_ab = y_intercept(lines(1,:),a,b);
    yz_bc = y_intercept(lines(2,:),b,c);
    yz_cd = y_intercept(lines(3,:),c,d);
    yz_da = y_intercept(lines(4,:),d,a);

    yz_ef = y_intercept(lines(5,:),e,f);
    yz_fg = y_intercept(lines(6,:),f,g);
    yz_gh = y_intercept(lines(7,:),g,h);
    yz_he = y_intercept(lines(8,:),h,e);

    y_ints = [yz_ab; yz_bc; yz_cd; yz_da; yz_ef; yz_fg; yz_gh; yz_he];


###########################################################################
function [int] = y_intercept(line,v1,v2);

    # find y = 0 for a line with params defined in lines
    # and end points given by v1 and v2

    m = line(1);
    b = line(2);

    # find x for y = 0;
    if m~=inf & m~=0
    x = -b/m;
    elseif m==inf
    x = b;  # recall if m=inf the x-val is in b
    elseif m==0
    x = inf;
    end


    # limits of line
    hx = max(v1(1), v2(1));
    lx  = min(v1(1), v2(1));

    if (x > hx) | (x < lx)
    int = [x,1];
    else;
    int = [x,0];
    end

###########################################################################
function lines = line_foot_print(a,b,c,d,e,f,g,h,i,j,k,l)
    # BEAM
    # calc equations for line a-b
    [lines(1,1),lines(1,2)] = line_param(a,b);
    # calc equations for line b-c
    [lines(2,1),lines(2,2)] = line_param(b,c);
    # calc equations for line c-d
    [lines(3,1),lines(3,2)] = line_param(c,d);
    # calc equations for line d-a
    [lines(4,1),lines(4,2)] = line_param(d,a);

    #DETECTOR
    # calc equations for line e-f
    [lines(5,1),lines(5,2)] = line_param(e,f);
    # calc equations for line f-g
    [lines(6,1),lines(6,2)] = line_param(f,g);
    # calc equations for line g-h
    [lines(7,1),lines(7,2)] = line_param(g,h);
    # calc equations for line h-e
    [lines(8,1),lines(8,2)] = line_param(h,e);

    #SAMPLE
    # calc equations for line i-j
    [lines(9,1),lines(9,2)] = line_param(i,j);
    # calc equations for line j-k
    [lines(10,1),lines(10,2)] = line_param(j,k);
    # calc equations for line k-l
    [lines(11,1),lines(11,2)] = line_param(k,l);
    # calc equations for line l-i
    [lines(12,1),lines(12,2)] = line_param(l,i);
    lines;
    #####################################################
