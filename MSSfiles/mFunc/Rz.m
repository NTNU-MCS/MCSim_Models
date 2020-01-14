function R = Rz(psi);
% R = Rz(psi) computes the Euler angle rotation matrix R in SO(3) for a 
% principal rotation of psi radians around the z-axis.
%
% Roger Skjetne - 30.01.2018

R = [cos(psi) -sin(psi) 0;
     sin(psi)  cos(psi) 0;
       0         0      1];