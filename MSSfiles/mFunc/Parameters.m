%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
%
% Parameters.m
%
% Project:	"Marine Cybernetics Simulator: MCSim"
%
% Abstract:	Simulation parameters, to be edited by user
% 
% Inputs:	-
%
% Outputs:	-
%	
% Calls:	-
%
% Author:	Øyvind Smogeli, Autumn 2002
%
% Changes:	Continuously by OS
%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 

% Choose which vessel to simulate
% Vessel 1: Supply vessel, Lpp 80m, B 17.4m, T 5.6 m
% Vessel 2: ? (Include your own vessels here)
% .....
Sim.vessno 				= 2; 			% 1: Supply , 2: CS3

% Global constants
global Glob
Glob.g 					= 9.81;			% Constant of gravity
Glob.rhow				= 1025;			% Density of seawater
Glob.p_0				= 101325;		% Atmospheric pressure
Glob.p_w				= 1230;			% Evaporation pressure of water
Glob.epsilon			= 1E-3;			% Small number

% Model reduction matrix
Glob.H36 = [1 0 0 0 0 0; 0 1 0 0 0 0 ; 0 0 0 0 0 1];

% Initial values for LF position, velocity and acceleration.
Sim.Eta_init			= [0 0 0 0 0 0];
Sim.Nu_init				= [0 0 0 0 0 0];
Sim.Nu_d_init			= [0 0 0 0 0 0];

%----------------------------------------------------------------------------------------------------------------------------------
% Simulation settings. 

% You can choose to run both LF and WF or either one alone. You may
% also choose which excitation forces to use. 
Sim.WF_model 			= 1;			% 1 for including WF model, 0 for not
Sim.LF_model 			= 1;			% 1 for including LF model, 0 for not
Sim.Prescribed_motion 	= 0;			% 1 for prescribing the vessel motion in Simulink (edit the file). ONLY FOR DEBUGGING!

% Load modules
Sim.Current_load		= 0;			% 1 for including current load model (= nonlinear damping), 0 for not
Sim.Wind_load			= 0;			% 1 for including wind load model, 0 for not
Sim.Wavedrift_load		= 1;			% 1 for including wavedrift load model, 0 for not. 
Sim.Force_plots			= 1;			% 1 for plotting the load time series, 0 for not. 

% DEBUGGING OPTIONS, all set to 1 for normal simulations
Sim.C_rb				= 1;			% 1 for including C_rb (rigid-body coriolis-centripetal forces), 0 for not
Sim.C_a					= 1;			% 1 for including C_a (added mass coriolis-centripetal forces), 0 for not
Sim.D_l					= 1;			% 1 for including D_l (linear damping forces), 0 for not
Sim.G					= 1;			% 1 for including G (restoring forces), 0 for not
Sim.Nu_r_C_a			= 1;			% 1 for using Nu_r in C_a, 0 for using Nu
Sim.Nu_r_D_l			= 1;			% 1 for using Nu_r in D_l, 0 for using Nu

% Reduced sampling rates for time-consuming routines, 1 second = 1 Hz should be adequate for most, but less for the thruster module
Sim.WF_sample_time		= 1;			% Sample time for the WF model in seconds
Sim.WD_sample_time		= 1;			% Sample time for the Wavedrift load model in seconds
Sim.currload_sample_time= 1;			% Sample time for the current load model in seconds
Sim.windload_sample_time= 1;			% Sample time for the wind load model in seconds
Sim.gust_sample_time	= 1;			% Wind gust model sample time in seconds
Sim.t_wave_sample_time	= 0.1;			% Thruster wave calculation sample time

% External force, for testing response, mostly for debugging
Sim.Ext_force 				= [0 0 0 0 0 0];	% Force and moment to be applied in vessel coordinate system
Sim.Ext_force_coord			= [0 0 0];			% Coordinates in vessel coordinate system at which force is applied
Sim.Ext_force_T				= 20;				% Time constant for ramping up external force

% Thruster simulation settings
Sim.Thrusters.Enable		= 0;		% 1 for including thruster dynamics in the simulation, 0 for not
										% If thrusters are disabled, the thrusters are represented by a simplified model.

% Thrust allocation enable
Sim.Thrust_alloc		= 0;			% 1 for enabling thrust allocation, 0 for disabling

% Sensor module enable
Sim.Ideal_meas			= 0;			% 1 for using ideal LF measurements, 0 for using the actual (LF + WF) measurements
Sim.Sensors				= 0;			% 1 for including sensor models, 0 for not. NB! Sensor module not fully implemented!

Sim.Sensor.Eta_noise	= 0;			% 3DOF noise (GPS + compass) for Eta, if Sim.Sensors == 1

% Observer enable
Sim.Vessel_Observer		= 2;			% 0 for disabling vessel observer and using measurements in feedback
										% 1 for enabling vessel observer but using measurements in feedback
										% 2 for enabling vessel observer and using estimations in feedback
										% NB! If set to zero or one, Sim.Sensors should be to to zero and Sim.Ideal_meas set to 1
% DP control system enable
Sim.Vessel_Control		= 0;			% 1 for enabling control system, 0 for disabling										
			
% Enabling/disabling of terms in controller. FOR DEBUGGING!
Sim.Vessel.HP.Controller.ref_ff   	= 1;       % 1 for enabling model reference feedforward term, 0 for disabling
Sim.Vessel.HP.Controller.int_act  	= 1;       % 1 for enabling integral action term, 0 for disabling
Sim.Vessel.HP.Controller.pd_act   	= 1;       % 1 for enabling PD-term, 0 for disabling

% Observer settings
Sim.Vessel.Observer.wave_filter   	= 1;       % 1 for enabling wave filtering in observer, 0 for disabling
Sim.Vessel.Observer.bias_time     	= 0;       % 1 for including bias time constant in observer, 0 for T = 0

%----------------------------------------------------------------------------------------------------------------------------------
% Settings for simplified thruster model, only for Sim.Thrusters == 0

Sim.Thrusters.thrust_sat	= 1;			% 1 for enabling saturation of thrust vector in simple thruster dynamics, 0 for not
Sim.Thrusters.thrust_time	= 1;			% 1 for enabling 1st order thruster dynamics, 0 for not

% Max thrust vector. Only used when Sim.Thrusters.thrust_sat == 1
Sim.Thrusters.Tau_max		= [10 10 12];	

% Time constants for 1st order thruster system, when Sim.Thrusters.thrust_time == 1
Sim.Thrusters.T_surge		= .5;		% s					
Sim.Thrusters.T_sway		= .7;		% s					
Sim.Thrusters.T_yaw			= .7;		% s					
										
%----------------------------------------------------------------------------------------------------------------------------------
% Settings for full thruster dynamics, only for Sim.Thrusters == 1

Sim.Thruster_enable		= [1 1 1 1 1];	% 1 for enabling, 0 for disabling. Only appliccable if Sim.Thrusters == 1;
Sim.nthrusters			= length(Sim.Thruster_enable);	% Number of thrusters. This must comply with the number of thruster units defined in the Simulink file. NB! ØK TIL 10

% Individual settings for each thruster unit
Sim.Thruster{1}.Wavecalc	= 1;		% Wave elevation and wave induced velocity calculations, 0 for excluding
Sim.Thruster{1}.Cc_loss		= 1;		% Cross coupling loss for ducted/open prop, transverse loss for tunnel thruster
Sim.Thruster{1}.Vent_loss	= 0;		% Ventilation and disc area losses
Sim.Thruster{1}.Wagner		= 0;		% "Wagner effect" during ventilation: Slower build-up than loss of thrust and torque
										% This is not fully functional yet, transient problems with changing rpm in a critical area

Sim.Thruster{1}.Hydro_damp	= 0;		% 1 for including hydrodynamic mass-damper dynamics in model, 0 for not. This is not very well tested....	
Sim.Thruster{1}.Damp_T_true = 0;		% 1 for using true thrust with losses in hydrodamic damper model, 0 for using nominal thrust
										
Sim.Thruster{1}.Fixed_rpm	= 0*60;		% Sets the thruster rpm fixed to this value if different from zero.
% This is ONLY for simulating thruster losses without thruster dynamics. NORMALLY = 0 !!!

Sim.Thruster{2} = Sim.Thruster{1};
Sim.Thruster{3} = Sim.Thruster{1};
Sim.Thruster{4} = Sim.Thruster{1};
Sim.Thruster{5} = Sim.Thruster{1};

% Sim.Thruster{2}.Hydro_damp	= 1;
% 
% Sim.Thruster{3}.Hydro_damp  = 1;
% Sim.Thruster{3}.Damp_T_true = 1;

%----------------------------------------------------------------------------------------------------------------------------------


%----------------------------------------------------------------------------------------------------------------------------------
% Set all control options and parameters in Controls.m										
%----------------------------------------------------------------------------------------------------------------------------------
	
%----------------------------------------------------------------------------------------------------------------------------------
% Environment settings

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

Par.Wave.time_change(1) = 365.16;          % Time to begin the change the wave
Par.Wave.time_change(2) = 365.53;          % Time to end the change the wave
for i=3:2:48
    Par.Wave.time_change(i)   = Par.Wave.time_change(i-2)+91.30;          % Time to begin the change the wave
    Par.Wave.time_change(i+1) = Par.Wave.time_change(i-1)+91.30;          % Time to end the change the wave
end

% Par.Wave.hs(1)			= 0.04;			% Significant wave heigth in sea state
% Par.Wave.omega_peak(1)	= 4.33;        	% Peak frequency component in sea state, 0 for expectation value based on Hs
% for i=2:25
%     Par.Wave.hs(i)		    = Par.Wave.hs(i-1)+0.0120833;			% Significant wave heigth in sea state 
%     Par.Wave.omega_peak(i)  = Par.Wave.omega_peak(i-1)-0.09125;         	% Peak frequency component in sea state, 0 for expectation value based on Hs
% end

Par.Wave.hs(1)			= 0.04;			% Significant wave heigth in sea state
Par.Wave.omega_peak(1)	= 4.33;        	% Peak frequency component in sea state, 0 for expectation value based on Hs
for i=2:25
    Par.Wave.hs(i)		    = Par.Wave.hs(i-1)+0.01708333333333;			% Significant wave heigth in sea state 
    Par.Wave.omega_peak(i)  = Par.Wave.omega_peak(i-1)-0.1388;         	% Peak frequency component in sea state, 0 for expectation value based on Hs
end

Par.Wave.psi_mean		= 180*pi/180;	% Mean direction of sea state (0 is aft, 180 stern)
Par.Wave.spread			= 4;			% Spreading factor for direction spectrum. 1 recommended by ITTC, 2 by ISSC. High value centers wave energy
Par.Wave.gamma			= 3.3;			% Gamma value for JONSWAP spectrum, if using this. Default is 3.3

Par.Wave.nwave  	    = 50;			% Number of waves <=nfreq*ndir
Par.Wave.nfreq  	    = 20;			% Number of frequency components in simulation
Par.Wave.ndir			= 10;			% Number of wave directions in simulation
Par.Wave.rand_freq		= 1;			% 1 for random frequencies
Par.Wave.rand_dir		= 1;			% 1 for random directions
Par.Wave.cutoff 		= 2.5;			% Cutoff frequency for spectrum = cutoff*omega_peak
Par.Wave.cutoff_dir		= 0;			% Cutoff direction for spectrum = cutoff_dir+-omega_peak

% Wind
Par.Wind.u10			= 0;			% 1 hour mean wind velocity at 10 meters. -1 for using the expectation value based on Hs and surface drag coefficient (also changes the direction)
Par.Wind.dir			= 30*pi/180;	% Mean wind direction (rad). Set equal to mean wave direction if vel = -1
Par.Wind.gust			= 0;			% 1 for including wind gust, 0 for not
Par.Wind.spectrum		= 2;			% 1 for Harris wind spectrum, 2 for NORSOK wind spectrum
Par.Wind.kappa			= 0.0026;		% Sea surface drag coefficient for calculation of mean wind speed. Harris wind spectrum: kappa = 0.0026.
Par.Wind.kappa_H		= 0.0026;		% Sea surface drag coefficient Harris wind spectrum.
Par.Wind.L				= 1800;			% Scaling length for Harris wind spectrum
Par.Wind.n 				= 0.468;		% NORSOK spectrum parameter
Par.Wind.nfreq			= 20;			% Number of frequency components in gust realization
Par.Wind.minfreq		= 1E-4;			% Minimum frequency in gust realization [Hz]
Par.Wind.maxfreq		= 1E-1;			% Maximum frequency in gust realization [Hz]
Par.Wind.var_dir		= 1;			% 1 for variation of direction, 0 for constant
Par.Wind.mju_dir		= 0;			% Gauss Markov process constant for varying direction. 0 for random walk.
Par.Wind.sat_dir		= 5*pi/180;		% Max variation in wind direction: direction = [dir - sat_dir, dir + sat_dir]
Par.Wind.wn_t_dir 		= 0.5;			% Direction white noise sample time
Par.Wind.wn_p_dir 		= 0.000005;		% Direction white noise power
Par.Wind.wn_s_dir 		= 6375;			% Direction white noise seed


% Current
% This is just for surface velocity. Define the current profile including depth below as Par.Current.Profile
% NB! DEFINE THE MEAN CURRENT DIRECTION AS FUNCTION OF TIME TO FACILITATE
% FAST CHANGES
Par.Current.vel 		= 0;			% Current velocity (m/s)
Par.Current.dir			= 160*pi/180;	% Current direction (rad)
Par.Current.var_vel		= 0;			% 1 for variation of velocity, 0 for constant
Par.Current.var_dir		= 0;			% 1 for variation of direction, 0 for constant
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

% Define current profile [xvel, yvel, depth]. If only using surface component, it is fully defined by the surface velocity and direction
Par.Current.Profile = ...
	[Par.Current.vel*[cos(Par.Current.dir), sin(Par.Current.dir)], 0;
	 Par.Current.vel*[cos(Par.Current.dir), sin(Par.Current.dir)], 5];

%----------------------------------------------------------------------------------------------------------------------------------
% Definitions
Par.Wave.spectypes{1} 	= 'ITTC';
Par.Wave.spectypes{2} 	= 'JONSWAP';
Par.Wave.spectypes{3} 	= 'Doubly Peaked';
Par.Wind.spectypes{1}	= 'Harrison';
Par.Wind.spectypes{2}	= 'NORSOK';

Sim.Sensor.Types{1}		= 'Off';
Sim.Sensor.Types{2}		= 'On';

Sim.Sensor.Feedback{1}	= 'WF + LF';
Sim.Sensor.Feedback{2}	= 'LF only';

Sim.VessObs.Types{1}	= 'Disabled, using measurement in feedback';
Sim.VessObs.Types{2}	= 'Enabled, but using measurement in feedback';
Sim.VessObs.Types{3}	= 'Enabled, using output in feedback';
