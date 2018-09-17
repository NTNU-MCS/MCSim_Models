% Readme.m
%
% Short introduction to MCSim release Beta 4
%
% CONTENTS
% - Introduction
% - Quickstart
% - Overview of the simulator
% - Instructions for use
% - Some pogramming rules/tips
%
% 
% THIS IS A BETA VERSION. EXPECT SOME PROBLEMS, AND HELP ME FIX THEM!!
%
% The Simulink code has been developed with Simulink 5.0, which is bundled with Matlab release 13.
%
% You need the gnc toolbox installed in your Matlab directory. Get this from
% http://www.marinecybernetics.com/software.htm
%
% NB!! SIMULINK BUG!
% From time to time Simulink will remove the links to the input ports of most of the lookup tables
% used in the simulator, producing a lot of error messages and bad results. If this happens, check
% the lookup tables in the simulator. The ones that have been troublesome so far, are:
% MCSim/Vessel_module/Vessel_dynamics/Vessel WF 6DOF/Wave component loop/6DOF transfer functions/TF tables
% MCSim/Vessel_module/Vessel_dynamics/Vessel LF 6DOF/Wavedrift load/Wavedrift force/Wave component loop/Force calculation/TF tables
% 
%
% 
% QUICKSTART
% - Set your working path to the MCSim root, where the file MCS.m is located.
% - Run the script MCS.m
% - Run the script Init.m
% - Open the Simulink file MCSim.mdl
% - Press play
% - Display results with the script Vessel_results.m
%
%
% OVERVIEW OF THE SIMULATOR
% 
% Environmental module
% 
% 	Wave model
% 	The sea state is specified by significant wave height Hs, mean wave direction, frequency spectrum,
% 	directional spectrum and the number of frequencies and direction in the spectral grid. ITTC, ISSC,
% 	Jonswap and Torsethaugen doubly peaked frequency spectra are included. To reduce the number of harmonic
% 	wave components, the components with little energy content may be discarded. This is an effective way
% 	of speeding up the simulation. A more realistic	sea state with few wave components can be achieved by
% 	choosing frequencies and directions randomly within each interval. The peak spectral frequency can be
% 	chosen manually, or an expectation value based on Hs (recommended by NORSOK) may be chosen. Plots of
% 	the wave spectrum and a realization of the sea surface are available as initialization options.
% 	
% 	Current model
% 	The surface current is defined by mean direction and magnitude. A stochastic variation of direction
% 	and magnitude may be chosen. No current profile is currently implemented, but a suggestion for structure
% 	is included.
% 	
% 	Wind model
% 	The mean wind direction and magnitude may be set manually, or expectation values based on Hs (recommended
% 	by NORSOK) may be chosen. Harrison or NORSOK wind gust spectra are added to the mean wind, and the wind
% 	direction may be chosen as a stochastic process.
% 
% Vessel module
%
% The vessel module is the core of the model, and contains several modules. The most important ones are described
% in short in the following.
% 
% 	Vessel Dynamics
% 
% 	The vessel dynamics consists of a low-frequency (LF) and a wave-frequency (WF) model, which may be enabled and
% 	disabled individually. If no actual dynamics is desired, a module for specifying a prescribed motion is included.
% 	
% 		LF module
% 		
% 		The LF module is based on the 6DOF kinematic equations as defined by Fossen. Asymptotic values for low-frequency
% 		added mass and linear damping are used. Nonlinear damping is included in the current load module. Wind load
% 		and Wavedrift load modules are also included.
% 				
% 		WF module
% 
% 		The WF module is currently based on 6DOF motion transfer functions. The harmonic motion due to each wave component
% 		is calculated and summed. This means that no external interaction force will affect the WF motion. A option for
% 		using load transfer functions will be added later, enabling true interaction between the vessel and external
% 		elements.
%
%		Prescribed motion
%	
%		If wanted, a vessel motion can be prescribed in this block and added to the LF and WF motion. This should only be
%		used for testing purposes.
% 				
% 	Thruster module
% 	
% 	Thruster dynamics may be enabled or disabled. In disabled mode, simplified thruster dynamics may be added in form of
%	saturation and 1st order systems. The desired thrust vector from the vessel control system is then used directly in
%	the vessel dynamics. If thrusters are enabled and the DP system is enabled, thrust allocation must also be enabled in
%	order to facilitate DP operation. 
% 	
% 		Thrust allocation
% 		
% 		Two thrust allocation options are available: Feedforward allocation (preset setpoint for each thruster) and
% 		thrust allocation for non-rotatable thrusters (generalized inverse). The second must be used with the DP system.
% 		
% 		Thruster units
% 		
% 		At the moment, 5 thruster units are available in the simulation. Each thruster unit has a mask with one parameter:
%		The thruster index k. This is a reference to the Par.Thrusters{k} parameter structure, which facilitates identical
%		code in each unit: If a change is made to the thruster unit, delete all the others,	replace them with copies of the
%		new one and update the index in the mask. Each unit may be enabled and disabled in the initialization.	A thruster
%		unit contains both local thruster control and thruster dynamics. The data structure Par.Thrusters is defined in
%		Thruster_config.m, which again refers to the thruster data defined in Thruster_data.m.
%
% 			Thruster control
% 			
% 			Four different local thruster control strategies are available:
% 			
% 			- Shaft speed PID control
% 			- Torque feedforward control
% 			- Power feedback control
% 			- Hybrid power/torque control
% 			
% 			Parameters for the controllers are set in Controls.m
%
% 			Thruster dynamics:
% 			
% 			Presently, the following nominal thruster models are available:
% 			
% 			- 4 quadrant model
% 			- 1 quadrant model
% 			- tunnel thruster
% 			
% 			The motor dynamics in the form of a 1st order system and shaft dynamics are included. Azimuth dynamics are not implemented
% 			properly. Thrust losses from to inline and transverse velocity fluctuations due to current, vessel motion and wave induced
% 			velocitites are calculated, as well as severe thrust losses due to ventilation and in-and-out-of water effects.
%
%			Simulation settings are set in Parameters.m 			
% 	
% 	Sensor module
% 	
% 	This is not fully implemented as of yet, but a choice of using only LF or LF and WF motion components in the feedback loop may be done.
%	A simple noise model for Eta is also included.
%
%	Vessel observer
%
% 	A nonlinear observer (with or without wave filtering) is included, and may be run in three modes:
% 	
% 	- Disabled
% 	- Enabled but outside the loop (without using the output as feedback to the DP controller)
% 	- Enabled and in the loop, giving feedback to the DP controller
% 	
% 	The wave filtering may be turned on and off, but the wave filtering frequency must be set manually. The observer requires quite a lot of
% 	tuning.
%
% 	Vessel controller
% 	
% 	The following vessel controllers are available:
%	
%	- Thrust vector feedforward (no feedback)
% 	- A simple PID DP controller
% 	- The same PID controller with online tuning options
% 	- Horizontal-plane DP controller with reference feedforward
% 	
% 	All parameters are set in Controls.m
% 	
%
% INSTRUCTIONS FOR USE
%
% To run a simulation, you need to do the following:
% 
% 1)	Set your working directory to the folder with this file (and the MCSim.mdl file).
%		You need to have the gnc toolbox path set in your Matlab search path.
%
% 2)	Run the script MCS.m and follow instructions to set the rest of the work path correctly.
%		To reset the Matlab path, run the script reset_path.m
%
% 3)	Set the simulation parameters in Parameters.m and control system parameters in Controls.m. 
%		Check carefully that you've set all options and parameters correctly for your simulation.
%		The thruster configuration is defined in Thruster_config.m, using the different thruster
%		units defined in Thruster_data.m
%
% 4)	Run the script Init.m to initialize the simulation. This reads all vessel and simulation
%		data into the workspace. If you want to, you could add this file as a callback function for
%		the MCSim model. This automatically initializes the simulation before it is run.
%
% 5)	Open the MCSim simulink file, and set the desired simulation settings (solver, time, step etc.).
%		Run the simulation. The simulator should work in both normal and accelerator mode. After the
%		simulation, relevant data is stored in a structure called Saved by the callback function Save_struct.
%		Edit the file to your own liking.
%
%		Check that no error messages have occured in the matlab command window. If one occurs, locate it
%		and try to fix it. There should normally be no error messages or warnings!
%
% 6)	Save your results with the script Savedata.m if wanted. The Saved structure created at simulation
%		end is then saved to a mat-file in the Results folder together with all the Simulation settings.
%		If you need to store other data or make different plots, just follow my setup or make one to your own liking.
%
% 7)	Plot results with Vessel_results.m, Thruster_results.m and Control_results.m. Edit the files as you like.	
%
% 8)	If you wish to plot an old data set, just load the mat-file (containg the Saved structure) to the workspace
%		and run Results.m
%
% DEFAULT SETTINGS
%
% The configuration files Thruster_config.m, Controls.m and Parameters.m included in the Beta4.0 distribution of MCSim
% are set up for a "typical" supply vessel equipped with three tunnel thrusters (two fore and one aft) and two ducted
% propulsors with fixed azimuth. Low-freqeuncy (LF) dynamics are simulated with loads from current, wind and wavedrift. The
% wave-frequency motion (WF) is added to the LF motion from motion transfer functions. A nonlinear observer with wave
% filtering provides feedback to a horizontal-plane DP controller with reference feedforward. The thrust allocation is a
% simple generalized inverse. No sensor models are enabled, so the feedback signals to the observer are noise free.
% Thruster dynamics with some loss effects are included, and hybrid torque/power controllers provide local thruster control.
% To run the simulation without thruster dynamics, set the flag Sim.Thrusters.Enable to 0. This will speed up the simulation.
%
% MAKING YOUR OWN APPLICATIONS
%
% There's basically one rule: KEEP THE MODULARITY. Try to make as few changes as possible to the existing code,
% and make your own applications modular. This will make debugging easier, and integration of many modules much
% simpler. If major changes are necessary to facilitate your application, come and see me first.
%
% Stay away from S-functions and Matlab-functions, as these slow down the simulation significantly (and makes
% it impossible to run in acceleration mode). S-functions written in C are allowed.
%
% If you need to delete part of the simulator for your own application, feel free to do so. Remember that most modules
% may be disabled by setting flags in the Parameters file, and then have no contribution to the CPU time.
%
% VESSEL DATA
% The current vessel data is for a supply vessel. If you want to specify other vessel data, you have to stick
% with the structures I have used. Make your OWN files, and specify a new vessel number. DON'T just overwrite
% the old data or try some other hack. And make sure you get the coordinate systems right.
%
% GOOD LUCK.
%
% Trondheim, 09/02 2004
% Oyvind Notland Smogeli
% smogeli@marin.ntnu.no