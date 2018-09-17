%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
% WaveParameters.m
%
% Project: ADPRC -	AMOS DP Research Cruise
%
% Abstract:	Wave parameters - to be edited by user
%
%
% Øyvind Smogeli
%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% Global constants
global Glob
Glob.g 					= 9.81;			% Constant of gravity
Glob.rho_w				= 1025;			% Density of seawater
Glob.p_0				= 101325;		% Atmospheric pressure
Glob.p_w				= 1230;			% Evaporation pressure of water
Glob.epsilon			= 1E-3;			% Small number


Sim.Wave.Enable = 1;   % 1 enabled, 0 disabled
Sim.Current.Enable = 1;  % 1 enabled, 0 disabled


% Waves
% Number of wave components = nfreq*ndir, but using the rand_dir and energylim options
% will reduce the number significantly. A large number of wave components slows the
% simulation down significantly, but for a realistic sea state you should use at least 100.
%Par.Wave.spectrum		= 3;			% Spectrum type: 1 = ITTC, 2 = JONSWAP, 3 = Doubly Peaked

%----------------------------------------------------------------------------------------------------------------------------------
% Environment settings

% Waves
% Number of wave components = nfreq*ndir, but using the rand_dir and energylim options
% will reduce the number significantly. A large number of wave components slows the
% simulation down significantly, but for a realistic sea state you should use at least 100.
Par.Wave.spectrum		= 2;			% Spectrum type: 1 = ITTC, 2 = JONSWAP, 3 = Doubly Peaked

Par.Wave.time_change(1) = 5000; %200;          % Time to begin the change the wave
Par.Wave.time_change(2) = 5000; %200.53;          % Time to end the change the wave
for i=3:2:48
    Par.Wave.time_change(i)   = Par.Wave.time_change(i-2)+439;  %91.30;          % Time to begin the change the wave
    Par.Wave.time_change(i+1) = Par.Wave.time_change(i-1)+439;   %91.30;          % Time to end the change the wave
end

% Par.Wave.hs(1)			= 0.04;			% Significant wave heigth in sea state
% Par.Wave.omega_peak(1)	= 4.33;        	% Peak frequency component in sea state, 0 for expectation value based on Hs
% for i=2:25
%     Par.Wave.hs(i)		    = Par.Wave.hs(i-1)+0.0120833;			% Significant wave heigth in sea state
%     Par.Wave.omega_peak(i)  = Par.Wave.omega_peak(i-1)-0.09125;         	% Peak frequency component in sea state, 0 for expectation value based on Hs
% end

Par.Wave.hs(1)			= 4;%6;%6;  %7.0;	%2;		% Significant wave heigth in sea state
Par.Wave.omega_peak(1)	= 0.6;%0.53;%.43;  %0.75;  %0.8;  % 0.5;  %.75;   % saet peak frequency very low when no waves     	% Peak frequency component in sea state, 0 for expectation value based on Hs
for i=2:25
    Par.Wave.hs(i)		    = Par.Wave.hs(i-1)+0.41668;%0.01708333333333;			% Significant wave heigth in sea state
    Par.Wave.omega_peak(i)  = Par.Wave.omega_peak(i-1)-0.010;         	% Peak frequency component in sea state, 0 for expectation value based on Hs
end

Par.Wave.psi_mean		= 150*pi/180;	% Mean direction of sea state  For 0 vessel heading 180 is head incident.
Par.Wave.spread			= 4;			% Spreading factor for direction spectrum. 1 recommended by ITTC, 2 by ISSC. High value centers wave energy
Par.Wave.gamma			= 3.3;			% Gamma value for JONSWAP spectrum, if using this. Default is 3.3

Par.Wave.nwave  	    = 10;			% Number of waves <=nfreq*ndir
Par.Wave.nfreq  	    = 5;			% Number of frequency components in simulation
Par.Wave.ndir			= 10;			% Number of wave directions in simulation
Par.Wave.rand_freq		= 1;			% 1 for random frequencies
Par.Wave.rand_dir		= 1;			% 1 for random directions
Par.Wave.cutoff 		= 2.5;			% Cutoff frequency for spectrum = cutoff*omega_peak
Par.Wave.cutoff_dir		= 0;			% Cutoff direction for spectrum = cutoff_dir+-omega_peak


% Current
% This is just for surface velocity. Define the current profile including depth below as Par.Current.Profile
% NB! DEFINE THE MEAN CURRENT DIRECTION AS FUNCTION OF TIME TO FACILITATE
% FAST CHANGES
Par.Current.vel 		= 0.5;			% Current velocity (m/s)
Par.Current.dir			= 180*pi/180;	% Current direction (rad)   for 0 vessel heading 180 is head incident.
Par.Current.var_vel		= 1;			% 1 for variation of velocity, 0 for constant
Par.Current.var_dir		= 1;			% 1 for variation of direction, 0 for constant
Par.Current.mju_vel		= 0;			% Gauss Markov process constant for varying velocity. 0 for random walk.
Par.Current.sat_vel		= 0.2;			% Max variation in current velocity: velocity = [vel - sat_vel, vel + sat_vel]
Par.Current.wn_t_vel 	= 0.5;			% Velocity white noise sample time
Par.Current.wn_p_vel 	= 0.000001;		% Velocity white noise power
Par.Current.wn_s_vel 	= 1234;			% Velocity white noise seed
Par.Current.mju_dir		= 0;			% Gauss Markov process constant for varying direction. 0 for random walk.
Par.Current.sat_dir		= 5*pi/180;		% Max variation in current direction: direction = [dir - sat_dir, dir + sat_dir]
Par.Current.wn_t_dir 	= 0.5;			% Direction white noise sample time
Par.Current.wn_p_dir 	= 0.000001;		% Direction white noise power
Par.Current.wn_s_dir 	= 4321;			% Direction white noise seed

Par.Current.vel2		= 0.05*0.5;			% Current velocity (m/s)
Par.Current.dir2			= 90*pi/180;	% Current direction (rad)   for 0 vessel heading 180 is head incident.

Current.Change.Time =5000; % time for a changing current (ad hoc solution for simulating a
% wave train ect.
TC_currentChange = 30; % time constant for the current change
% Define current profile [xvel, yvel, depth]. If only using surface component, it is fully defined by the surface velocity and direction
Par.Current.Profile = ...
    [Par.Current.vel*[cos(Par.Current.dir), sin(Par.Current.dir)], 0;
    Par.Current.vel*[cos(Par.Current.dir), sin(Par.Current.dir)], 5];

Par.Current.Profile2 = ...
    [Par.Current.vel2*[cos(Par.Current.dir2), sin(Par.Current.dir2)], 0;
    Par.Current.vel2*[cos(Par.Current.dir2), sin(Par.Current.dir2)], 5];



%----------------------------------------------------------------------------------------------------------------------------------
% Definitions
Par.Wave.spectypes{1} 	= 'ITTC';
Par.Wave.spectypes{2} 	= 'JONSWAP';
Par.Wave.spectypes{3} 	= 'Doubly Peaked';
Par.Wind.spectypes{1}	= 'Harrison';
Par.Wind.spectypes{2}	= 'NORSOK';
