% Init Gunnerus
% loads all necessary files for running Gunnerus on DP
% This program uses the MSS toolbox DP Motion RAO template

% ________________________________________________________________
%
% MSS HYDRO/GNC is a Matlab toolbox for guidance, navigation and control.
% The toolbox is part of the Marine Systems Simulator (MSS).
%
% Copyright (C) 2008 Thor I. Fossen and Tristan Perez
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% 
% E-mail: contact@marinecontrol.org
% URL:    <http://www.marinecontrol.org>

clear all
close all

% sample time [s]
Ts = 0.01;
% length of simulation [s]
Tsim = 200;  


% Load data from ShipX/VeRes (You can load data for another ship as well, if you have that...)
% LoadGummerus;
load GunnerusMan_input.mat
load GunnerusMan_inputABC.mat

% Load environmental parameters
EnvParameters;

ki=0;
kp=0;

k1=2*kp;
k2=7*kp;

Mb=1;
Mb_hat=2;

gamma=100;