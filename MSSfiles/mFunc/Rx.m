function R = Rx(phi);
% R = Rx(phi) computes the Euler angle rotation matrix R in SO(3) for a 
% principal rotation of phi radians around the x-axis.
%
% Roger Skjetne - 30.01.2018

R = [ 1     0          0     ; 
      0   cos(psi)  -sin(psi);
      0   sin(psi)   cos(psi)];
