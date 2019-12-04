function W = IMUsensormtrx(r)
% W = IMUsensormtrx(r) computes the 3x12 IMU configuration matrix 
% that relates the IMU sum of linear acceleration, angular acceleration, 
% and angular rate cross products to 3x1 measurement vector, that is, 
% xi = a_lin + S(r)'*alpha + H(r)*omega_bar where a_lin is the linear
% acceleration, alpha is the angular acceleration, and omega_bar is a 6x1
% vector of angular rate cross product terms.
% 
% For more info, see:
%  Kjerstad, Ø. K. and R. Skjetne, “Disturbance rejection by acceleration 
%  feedforward for marine surface vessels.” IEEE Access, Vol.4, 
%  pp.2656-2669, 2016. 
%
% Author:   Roger Skjetne
% Date:     14th November 2017
% Revisions: 
% ________________________________________________________________
%


S = [    0  -r(3)   r(2)
      r(3)     0   -r(1)
     -r(2)   r(1)     0 ];
 
H = [   0  -r(1) -r(1) r(2) r(3)   0
     -r(2)    0  -r(2) r(1)   0  r(3)
     -r(3) -r(3)    0    0  r(1) r(2)];

W = [eye(3) S' H];
