% Idealized 3DOF maneuvering model of a ship in currents.
%
% - Steady uniform current
% - Diagonal mass and added mass matrix
% - Linear + quadratic diagonal damping matrices
%
% function vessel_model() located at end of file. 
% 
%    Copyright (C) 2020: 	NTNU, Trondheim, Norway
%    Licensed under GPL-3.0-or-later
%    Created:  	<6-Nov-2020>	<Mathias Marley>
%    Revised:	<date>	<author> <description>
%

%% Preliminaries
clc; clear all; close all;
set(0,'defaultAxesFontSize',12)
set(0,'defaultFigureColor','w');
set(0,'defaultLineLineWidth',1.5);
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
plotx = 1000; ploty=600;

%% Model parameters
% Input parameters
L = 100; %Ship length
u0 = 6; v0 = 2; r0 = 2*pi/180; %reference velocities for linearized damping

% Calculated mass parameters
B = L/10; D = L/20; %ship length, breadth, draft [m]. 
m = L*B*D*1000; %rigid body mass [kg]
lr = L/4; %inertia radius in yaw
Mrb = diag([m m m*lr^2]); %Rigid body mass matrix
a11 = 0.05*m; %added mass surge [kg] 
a22 = 0.3*m; %added mass sway [kg] 
a33 = a22*L^2/4; %total yaw inertia (rigid body + added mass)
Ma = diag([a11 a22 a33]); %Added mass matrix

% Calculated damping coefficients
Cdx = 0.5; Cdy = 1; %quadratic damping coefficients
d11 = 0.5*B*D*Cdx*u0*1000; %linear damping
d11q = d11/u0; %quadratic damping
d22 = 0.5*L*D*Cdy*v0*1000;
d22q = d22/v0;
d33 = r0*Cdy*D*L^4/64*1000;
d33q = d33/r0;
Dl = diag([d11 d22 d33]); %Linear damping matrix

% Store model parameters in structure
param.Ma = Ma; param.Mrb = Mrb; param.m=m;
param.a11 = a11; param.a22 = a22;
param.Dl = Dl; param.d11q=d11q; param.d22q=d22q; param.d33q = d33q;

% Critical velocity for directional stability:
ucrit = sqrt(d22*d33/(m+a11)/(a22-a11));

%% Simulation example
% No controller. Constant actuator force in surge, step function in yaw. 
% Forward speed equal to critical speed for directional stability. 

Uc = 0.1; %[m/s], current velocity
betac = pi/4; %[rad], current direction in global reference frame

dt = 1; tvec = 0:dt:1800; 
ud = tvec*0+ucrit*1; %desired surge velocity

xvec = tvec*0; yvec=tvec*0; %ship position
uvec = ud; vvec = tvec*0; %surge, sway velocity
psivec = tvec*0; rvec = tvec*0; %yaw and yaw rate

Fuvec = ud.^2*d11q+ud.*d11; %actuator force in surge (constant)
Fvvec = 0*tvec; %no sway force

%actuator force in yaw (step function)
Frvec = tvec*0; 
Frvec(tvec>100)=d33*pi/180;
Frvec(tvec>700)=-d33*pi/180;
Frvec(tvec>1300)=0;

%% Time simulation
for i=1:length(tvec)-1
p = [xvec(i); yvec(i); psivec(i)];
nu = [uvec(i); vvec(i); rvec(i)];
F = [Fuvec(i); Fvvec(i); Frvec(i)];

%Forward euler time integration
[dp, dnu] = vessel_model(p,nu,F,Uc,betac,param);    
p = p+dp*dt; 
p(3)=atan2(sin(p(3)),cos(p(3))); %map to [-pi pi]
nu = nu+dnu*dt;

%Store results
xvec(i+1)=p(1); yvec(i+1)=p(2); psivec(i+1)=p(3); 
uvec(i+1)=nu(1); vvec(i+1)=nu(2); rvec(i+1)=nu(3);
end

%% Post-process
Uvec = sqrt(uvec.^2+vvec.^2); %Total speed
ucvec = Uc*cos(betac-psivec); %Body-fixed current velocities
vcvec = Uc*sin(betac-psivec); %Body-fixed current velocities
urvec = uvec-ucvec; vrvec=vvec-vcvec; %Relative velocities
betavec = atan2(vvec,uvec); %Drift angle
betarvec = atan2(vrvec,urvec); %Fluid relative drift angle
xpvec = -vvec./rvec; %pivot point (as defined)
xprvec = -vrvec./rvec; %pivot point defined versus relative velocity
xpest = (m+a11)*uvec/d22; % Estimate of steady-state pivot point from
%surge velocity (assuming linear damping dominates). 

%% Plots
close all
figure('Position',[300 100 plotx ploty])
plot(yvec,xvec); axis equal; hold on
xlabel('$y$ [m]'); ylabel('$x$ [m]');
title('Trajectory')

figure('Position',[300 100 plotx ploty])
subplot(2,1,1)
plot(tvec,uvec); hold on;
plot(tvec,urvec);
legend('$u$','$u_r$')
xlabel('Time [s]'); ylabel('[m/s]');
title('Surge velocity')
subplot(2,1,2)
plot(tvec,vvec); hold on
plot(tvec,vrvec)
legend('$v$','$v_r$')
xlabel('Time [s]'); ylabel('[m/s]');
title('Sway velocity')

figure('Position',[300 100 plotx ploty])
subplot(2,1,1)
plot(tvec,psivec*180/pi); hold on;
xlabel('Time [s]'); ylabel('$\psi$ [deg]');
title('Yaw angle')
subplot(2,1,2)
plot(tvec,rvec*180/pi); hold on
xlabel('Time [s]'); ylabel('$r$ [deg/s]');
title('Yaw velocity')

figure('Position',[300 100 plotx ploty])
subplot(2,1,1)
plot(tvec,Uvec); hold on;
xlabel('Time [s]'); ylabel('$U$ [m/s]');
title('Total speed')
ylim([0 ceil(max(Uvec))])
subplot(2,1,2)
plot(tvec,betavec*180/pi); hold on;
plot(tvec,betarvec*180/pi); hold on;
xlabel('Time [s]'); ylabel('[deg]');
legend('$\beta$','$\beta_r$')
title('Drift angle')

figure('Position',[300 100 plotx ploty])
subplot(2,1,1)
plot(tvec,xpvec); hold on;
plot(tvec,xprvec); hold on;
plot(tvec,xpest+tvec*0);
xlabel('Time [s]'); ylabel('[m]');
title('Pivot point')
legend('$x_p=-v/r$','$x_{pr}-v_r/r$','Estimate of steady-state')
ylim([0 L])

subplot(2,1,2)
plot(tvec,rvec*180/pi); hold on
xlabel('Time [s]'); ylabel('$r$ [deg/s]');
title('Yaw velocity')



%%

function [dp, dnu] = vessel_model(p,nu,F,Uc,betac,param)
% Maneuvering model of simplified ship. Assumes diagonal rigid body mass, 
% added mass and damping matrices (centripetal forces not included).
% Diagonal linear + quadratic damping. 
%
% Inputs:
%    p - [x y psi]' position and heading vector [m m rad]
%    nu - [u v r]' surge, sway and yaw velocity vector [m m rad/s]
%    F - [Fu Fv Fr] external force vector (e.g. actuator forces) [N N Nm]
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
%    Created:  	<06-Nov-2020>	<Mathias Marley>
%    Revised:	<date>	<author> <description>
%    		<date>	<author> <description>

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




