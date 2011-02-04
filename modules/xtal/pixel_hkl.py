"""
Pixel HKL

Authors/Modifications:
-----------------------

Notes:
------

"""
###########################################################

import numpy as num

###########################################################
""""
lam = 1.2;
k = 2*pi/lam;

# sample det distance (mm)
y = 5000 ;

# pixel dims (mm)
del_pixel = 0.175;

# pixel coords
dx = 500;
dz = 100;

# vector pointing to pixel
p = [dx*del_pixel, y, dz*del_pixel];

# normalized
pn = p/sqrt((dx*del_pixel)^2 + y^2 + (dz*del_pixel)^2);
v_mag(pn);


# get rot matricies for detector arm
nu = 0
del = 0

P = [  1,       0,         0    ;
       0,   cosd(nu), -sind(nu) ;
       0,   sind(nu),  cosd(nu) ];
    
D = [  cosd(del),   sind(del),   0 ;
      -sind(del),   cosd(del),   0 ;
          0     ,       0    ,   1 ];
       

# calc kr for center of det
kr = P*D*[0;1;0];
kr = k*kr
#should be same as
#kr = k*[ sind(del)           ; 
#         cosd(nu)*cosd(del)  ; 
#         sind(nu)*cosd(del) ];

# calc Q for center
ki = k*[ 0 ;
         1 ;
         0 ];      
Qm = kr - ki
tth = v_angle(ki,kr, [1 1 1 90 90 90])

# calc kr or rotated pixel vector 
kr_p = P*D*pn;
kr_p = k*kr_p

# calc Q for pixel
Qm_p = kr_p - ki
tth_p = v_angle(ki,kr_p, [1 1 1 90 90 90])

# now we can calc hkl for the center and the pixel 
# given Z and UB
# [Z,Qm,ki,kr,tth] = calc_Z_and_Q(angles,k);
#[U, UB] = calc_U(h_1, lam_1, angle_1, h_2, lam_2, angle_2, cell);
#
# h = (1/(2*pi))*inv(UB)*inv(Z)*Qm_p

########################################################################
########################################################################
#HKL scan for goethite (10L) rod, specfile: gt_6_april07_b(fake).spc, scan = 1 0 2.51379
#############################
# #S 2  hklscan  0.999939 0.999939  -8.56988e-05 -8.56988e-05  2.06 3.94  29 0.5
# #D Mon Apr 23 17:06:40 2007
# #T 0.5  (Seconds)
# #G0 0 0 1 -0.0001806167054 0.0002695003495 1 0 0 0 0 0 0 50.031 0 0 1 4 2 5 4 0
# #G1 3.023 4.609 9.963 90 90 90 2.078460241 1.363242636 0.6306519429 90 90 90 0 0 6 1 -1 3 -0.0002 -0.008 -2.59 1.9435 36.2778 15.309 26.6488 13.3275 9.021 -26.4045 13.88 -0.10275 1.0339 1.0339
# #G3 2.070853321 -0.09753606408 0.02949524505 0.1533409954 1.357964195 -0.03015295668 -0.08972198422 0.06964512884 0.6292397817
# #G4 0.9999385052 -8.569884087e-05 1.940140169 1.0339 4.000241947 7.548864431 -1.522131988 22.89000882 87.91053344 59.54799609 71.72405792 12.00256456 3.7875 45.7 4 0 0 0 12 0 0 0 0.2 122.798 0 0 0 0 -180 -180 -180 -180 -180 -180 -180 -180 -180 0 0 0 0 0 0 0 0 0 0 0 0
# #Q 0.999939 -8.56988e-05 1.94014
# #P0 21.675099 10.833 10.34 4.272 7.5422504 -1.2597501 15.1823 13.5044
# #P1 8.6210004 0 0.0001072 0.064599998 8.1185004 0 499.93943 0
# #P2 499.93943 -1816.025 179.96889 0 199.84443 0 2000.0779 0
# #P3 2000.0779 0 2000.0779 0 2000.0779 0.065399998 -0.1552 -0.1538
# #P4 0 0 0 0 0 0 0 -5655.4269
# #P5 0 
# #ATTEN 0 0 0 0
# #ALP_BET 4.00024 7.54886
# #ENERGY 11991.9
# #N 15
# #L H  K  L  del  eta  chi  phi  nu  mu  Alpha  Beta  Epoch  io  i1  Bicron  AmpTek_sc  ROI1  ROI2  ROI3  Seconds  IROI
# 0.999939 -8.56988e-05 2.06    21.8112    10.9035    10.3590     4.4855     8.3050    -1.2860  3.9997   8.2706  510 40 22 0 0 0 0 0 0.509996 0
# 0.999939 -8.56988e-05 2.12483    21.8852    10.9385    10.3700     4.6160     8.7168    -1.3003  3.9993   8.6601  512 40 23 0 0 0 0 0 0.509997 0
# 0.999939 -8.56988e-05 2.18966    21.9572    10.9755    10.3810     4.7540     9.1298    -1.3138  4.0002   9.0487  515 43 25 0 0 0 0 0 0.509996 0
# 0.999939 -8.56988e-05 2.25448    22.0292    11.0105    10.3940     4.9010     9.5435    -1.3308   3.998   9.4412  517 48 29 0 0 0 0 0 0.509996 0
# 0.999939 -8.56988e-05 2.31931    22.1011    11.0530    10.4060     5.0545     9.9590    -1.3453  3.9997   9.8307  519 48 30 0 0 0 0 0 0.509996 0
# 0.999939 -8.56988e-05 2.38414    22.1683    11.0820    10.4190     5.2180    10.3745    -1.3613  3.9977   10.223  521 42 23 0 0 0 0 0 0.509996 0
# 0.999939 -8.56988e-05 2.44897    22.2413    11.1175    10.4335     5.3880    10.7918    -1.3748  3.9999   10.615  523 48 29 0 0 0 0 0 0.509996 0
# 0.999939 -8.56988e-05 2.51379    22.3113    11.1500    10.4470     5.5685    11.2093    -1.3910   3.999   11.008  525 48 30 0 0 0 0 0 0.509996 0
# 0.999939 -8.56988e-05 2.57862    22.3787    11.1850    10.4620     5.7560    11.6275    -1.4060  4.0002   11.399  527 47 30 0 0 0 0 0 0.509996 0
# 0.999939 -8.56988e-05 2.64345    22.4432    11.2180    10.4770     5.9515    12.0480    -1.4228  3.9995   11.794  529 48 30 0 0 0 0 0 0.509997 0
###########################################################
# Give the hkl of interest
# Do the reverse converions from h to Pixels
h = [ 1.0021;-0.0028;2.5374]

# arrays from spec file
Uarray = [3.023 4.609 9.963 90 90 90 2.078460241 1.363242636 0.6306519429 90 90 90 0 0 6 1 -1 3 -0.0002 -0.008 -2.59 1.9435 36.2778 15.309 26.6488 13.3275 9.021 -26.4045 13.88 -0.10275 1.0339 1.0339];
angles = [22.3113    11.1500    10.4470     5.5685    11.2093    -1.3910 ];
# get the angles needed for the geometry calcs
# first parse the data needed to calc the UB matrix.
cell = [Uarray(1:6)];     # cell parameters (real space)
h_1 = [Uarray(13:15)];    # hkl of OR1
h_2 = [Uarray(16:18)];    # hkl of OR2 

# note see calc_U for how these should be parsed
angle_1 = [Uarray(19:24)];   # motor settings for OR1
angle_2 = [Uarray(25:30)];   # motor settings for OR2
lam_1 = Uarray(31);          # lambda for OR1  
lam_2 = Uarray(32);          # lambda for OR2 

lam = 1.0339;
k = 2*pi/lam;
# sample det distance (mm)
y =1500 ;
# pixel dims (mm)
del_pixel = 0.175;
#center pixels 
CenPix_x = 235;
CenPix_y = 85;
# now we can calc hkl for the center and the pixel 
# given Z and UB
 [Z,Qm,ki,kr,tth] = calc_Z_and_Q(angles,k);
[U, UB] = calc_U(h_1, lam_1, angle_1, h_2, lam_2, angle_2, cell);

#matrix of conversion of Pixel to h
M_PQ = (1/(2*pi))*inv(UB)*inv(Z);
Qm_p = inv(M_PQ)*h;
kr_p = Qm_p + ki;
kr_p = kr_p/k;
pn = inv(D)*inv(P)*kr_p;
mod_p = y/pn(2);
p = pn*mod_p;
dx = p(1)/del_pixel
dz = p(3)/del_pixel
Pix_x = CenPix_x + dz
Pix_y = CenPix_y + dx

########################################################################
########################################################################
#HKL scan for goethite (10L) rod, specfile: gt_6_april07_b(fake).spc, scan = 1 0 2.51379
#############################
# #S 2  hklscan  0.999939 0.999939  -8.56988e-05 -8.56988e-05  2.06 3.94  29 0.5
# #D Mon Apr 23 17:06:40 2007
# #T 0.5  (Seconds)
# #G0 0 0 1 -0.0001806167054 0.0002695003495 1 0 0 0 0 0 0 50.031 0 0 1 4 2 5 4 0
# #G1 3.023 4.609 9.963 90 90 90 2.078460241 1.363242636 0.6306519429 90 90 90 0 0 6 1 -1 3 -0.0002 -0.008 -2.59 1.9435 36.2778 15.309 26.6488 13.3275 9.021 -26.4045 13.88 -0.10275 1.0339 1.0339
# #G3 2.070853321 -0.09753606408 0.02949524505 0.1533409954 1.357964195 -0.03015295668 -0.08972198422 0.06964512884 0.6292397817
# #G4 0.9999385052 -8.569884087e-05 1.940140169 1.0339 4.000241947 7.548864431 -1.522131988 22.89000882 87.91053344 59.54799609 71.72405792 12.00256456 3.7875 45.7 4 0 0 0 12 0 0 0 0.2 122.798 0 0 0 0 -180 -180 -180 -180 -180 -180 -180 -180 -180 0 0 0 0 0 0 0 0 0 0 0 0
# #Q 0.999939 -8.56988e-05 1.94014
# #P0 21.675099 10.833 10.34 4.272 7.5422504 -1.2597501 15.1823 13.5044
# #P1 8.6210004 0 0.0001072 0.064599998 8.1185004 0 499.93943 0
# #P2 499.93943 -1816.025 179.96889 0 199.84443 0 2000.0779 0
# #P3 2000.0779 0 2000.0779 0 2000.0779 0.065399998 -0.1552 -0.1538
# #P4 0 0 0 0 0 0 0 -5655.4269
# #P5 0 
# #ATTEN 0 0 0 0
# #ALP_BET 4.00024 7.54886
# #ENERGY 11991.9
# #N 15
# #L H  K  L  del  eta  chi  phi  nu  mu  Alpha  Beta  Epoch  io  i1  Bicron  AmpTek_sc  ROI1  ROI2  ROI3  Seconds  IROI
# 0.999939 -8.56988e-05 2.06    21.8112    10.9035    10.3590     4.4855     8.3050    -1.2860  3.9997   8.2706  510 40 22 0 0 0 0 0 0.509996 0
# 0.999939 -8.56988e-05 2.12483    21.8852    10.9385    10.3700     4.6160     8.7168    -1.3003  3.9993   8.6601  512 40 23 0 0 0 0 0 0.509997 0
# 0.999939 -8.56988e-05 2.18966    21.9572    10.9755    10.3810     4.7540     9.1298    -1.3138  4.0002   9.0487  515 43 25 0 0 0 0 0 0.509996 0
# 0.999939 -8.56988e-05 2.25448    22.0292    11.0105    10.3940     4.9010     9.5435    -1.3308   3.998   9.4412  517 48 29 0 0 0 0 0 0.509996 0
# 0.999939 -8.56988e-05 2.31931    22.1011    11.0530    10.4060     5.0545     9.9590    -1.3453  3.9997   9.8307  519 48 30 0 0 0 0 0 0.509996 0
# 0.999939 -8.56988e-05 2.38414    22.1683    11.0820    10.4190     5.2180    10.3745    -1.3613  3.9977   10.223  521 42 23 0 0 0 0 0 0.509996 0
# 0.999939 -8.56988e-05 2.44897    22.2413    11.1175    10.4335     5.3880    10.7918    -1.3748  3.9999   10.615  523 48 29 0 0 0 0 0 0.509996 0
# 0.999939 -8.56988e-05 2.51379    22.3113    11.1500    10.4470     5.5685    11.2093    -1.3910   3.999   11.008  525 48 30 0 0 0 0 0 0.509996 0
# 0.999939 -8.56988e-05 2.57862    22.3787    11.1850    10.4620     5.7560    11.6275    -1.4060  4.0002   11.399  527 47 30 0 0 0 0 0 0.509996 0
# 0.999939 -8.56988e-05 2.64345    22.4432    11.2180    10.4770     5.9515    12.0480    -1.4228  3.9995   11.794  529 48 30 0 0 0 0 0 0.509997 0
###########################################################
# convert pixel to hkl from the center pixel
#pixels of interest
Pix_x = 255
Pix_y = 95

#########
Uarray = [3.023 4.609 9.963 90 90 90 2.078460241 1.363242636 0.6306519429 90 90 90 0 0 6 1 -1 3 -0.0002 -0.008 -2.59 1.9435 36.2778 15.309 26.6488 13.3275 9.021 -26.4045 13.88 -0.10275 1.0339 1.0339];
angles = [22.3113    11.1500    10.4470     5.5685    11.2093    -1.3910 ];
# get the angles needed for the geometry calcs
# first parse the data needed to calc the UB matrix.
cell = [Uarray(1:6)];     # cell parameters (real space)
h_1 = [Uarray(13:15)];    # hkl of OR1
h_2 = [Uarray(16:18)];    # hkl of OR2 

# note see calc_U for how these should be parsed
angle_1 = [Uarray(19:24)];   # motor settings for OR1
angle_2 = [Uarray(25:30)];   # motor settings for OR2
lam_1 = Uarray(31);          # lambda for OR1  
lam_2 = Uarray(32);          # lambda for OR2 

lam = 1.0339;
k = 2*pi/lam;

# sample det distance (mm)
y =1500 ;

# pixel dims (mm)
del_pixel = 0.175;

#center pixels 
CenPix_x = 235;
CenPix_y = 85;
# pixel coords
dx = (195 - CenPix_y) - (195 - Pix_y)
dz = Pix_x  - CenPix_x

# vector pointing to pixel
p = [dx*del_pixel; y; dz*del_pixel];

# normalized
pn = p/sqrt((dx*del_pixel)^2 + y^2 + (dz*del_pixel)^2);
v_mag(pn);

# get rot matricies for detector arm
nu = angles(5);
del = angles(1);

P = [  1,       0,         0    ;
       0,   cosd(nu), -sind(nu) ;
       0,   sind(nu),  cosd(nu) ];
    
D = [  cosd(del),   sind(del),   0 ;
      -sind(del),   cosd(del),   0 ;
          0     ,       0    ,   1 ];
       

# calc kr for center of det
kr = P*D*[0;1;0];
kr = k*kr;
#should be same as
#kr = k*[ sind(del)           ; 
#         cosd(nu)*cosd(del)  ; 
#         sind(nu)*cosd(del) ];

# calc Q for center
ki = k*[ 0 ;
         1 ;
         0 ];      
Qm = kr - ki;
tth = v_angle(ki,kr, [1 1 1 90 90 90]);

# calc kr or rotated pixel vector 
kr_p = P*D*pn;
kr_p = k*kr_p;

# calc Q for pixel
Qm_p = kr_p - ki;
tth_p = v_angle(ki,kr_p, [1 1 1 90 90 90]);

# now we can calc hkl for the center and the pixel 
# given Z and UB
 [Z,Qm,ki,kr,tth] = calc_Z_and_Q(angles,k);
[U, UB] = calc_U(h_1, lam_1, angle_1, h_2, lam_2, angle_2, cell);
h = (1/(2*pi))*inv(UB)*inv(Z)*Qm_p

"""


