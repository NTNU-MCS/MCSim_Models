function R = Ry(theta);
% R = Ry(theta) computes the Euler angle rotation matrix R in SO(3) for a 
% principal rotation of theta radians around the y-axis.
%
% Roger Skjetne - 30.01.2018

R = [ cos(psi)  0  sin(psi);
        0       1   0      ;
     -sin(psi)  0  cos(psi)];
       