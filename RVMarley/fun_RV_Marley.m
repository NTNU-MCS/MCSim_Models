function [dp, dnu] = fun_RV_Marley(p,nu,F,Uc,betac,param)
% Maneuvering model of an idealized ship named RV Marley (RV=Research Vessel)
% Assumes diagonal rigid body mass, added mass and damping matrices. 
% Centripetal forces not included in Coriolis matrix, so will give
% incorrect result without warning if non-diagonal mass matrices are used. 
% Diagonal linear + quadratic damping. 
%
% Inputs:
%    p - [x y psi]' position and heading vector [m m rad]
%    nu - [u v r]' surge, sway and yaw velocity vector [m m rad/s]
%    F - [Fu Fv Fr]' external force vector (e.g. actuator forces) [N N Nm]
%    Uc - current speed [m/s]
%    betac - current direction [rad]
%    param - structure containing ship parameters
%       param.Ma - diagonal added mass matrix 
%       param.Mrb - diagonal rigid body mass matrix 
%       param.Dl - diagonal linear damping matric 
%       param.d11q, param.d22q, param.d33q - quadratic damping parameters
%    
% Outputs:
%    dp - time derivative of p vector
%    dnu - time derivative of nu vector
%
%    Copyright (C) 2020: 	NTNU, Trondheim, Norway
%    Licensed under GPL-3.0-or-later
%    Created:  	11-Nov-2020	Mathias Marley
%    Revised:  	<date>	<author> <description>

%Unwrap surge, sway, yaw velocities
u=nu(1); v=nu(2); r = nu(3); psi = p(3);

%Calculate current and relative velocities:
uc = Uc*cos(betac-psi); ur = u-uc;
vc = Uc*sin(betac-psi); vr = v-vc;
nu_c = [uc vc 0]'; %Current velocity vector
dnu_c = Uc*[sin(betac-psi); cos(betac-psi); 0]*r; %current derivative
nu_r = nu-nu_c; %relative velocity

%Unrwap mass parameters
Mrb = param.Mrb; Ma=param.Ma; M = Mrb+Ma;
m = Mrb(1,1); a11 = Ma(1,1); a22=Ma(2,2);

%Calculate Coriolis forces
Fcoriolis_rb = [-m*v*r;m*u*r;0]; %rigid body
Fcoriolis_a = [-a22*vr*r; a11*ur*r; (a22-a11)*ur*vr]; %added mass

%Calculate damping matrix
d11q = param.d11q; d22q=param.d22q; d33q = param.d33q;
Dnl = diag([d11q*abs(ur);d22q*abs(vr);d33q*abs(r)]);
D = param.Dl+Dnl;

%Calculate time derivatives
dnu = inv(M)*(F-D*nu_r-Fcoriolis_rb-Fcoriolis_a+Ma*dnu_c);
dp = [cos(psi) -sin(psi) 0; sin(psi) cos(psi) 0; 0 0 1]*nu;

end




