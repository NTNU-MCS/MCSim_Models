function Par_Thruster = Thruster_config
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
% Par_Thruster = Thruster_config
%
% Project:	"Marine Cybernetics Simulator: MCSim"
%
% Abstract:	Loads thruster data	for simulation, defines the thruster positions and
%			orientations. As many thrusters as defined in the Simulink MCSim file
%			must be initialized. If not all are used, disable them in the Parameters.m
%			file.
% 
% Inputs:	- 
%
% Outputs:	Par.Thruster
%
% Calls:	Thruster_data
%						
% Author:	Oyvind Smogeli, Summer 2003
%
% Revision:	-
%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 

% Load all thruster defintions into structure "Thrusters"
global Glob
Thruster_data;

% Define the thrusters
% Thruster types:
% 1 - Ducted main prop, 4 quadrant data
% 2 - Unducted main prop, 4 quadrant data
% 3 - Unducted main prop, 1-quadrant data
% 4 - Unducted main prop, nominal data (only KT0 and KQ0, no Va dependence)
% 5 - Tunnel thruster, bow mounted
% 6 - CPP propeller, no pitch control implemented yet

% Thruster 1
Par_Thruster{1} 		= Thrusters{1};		% Ducted main prop
Par_Thruster{1}.r_p		= [-37.8 -5 4]';	% thruster position(x,y,z)
Par_Thruster{1}.alpha	= 0;				% thruster orientation

% Thruster 2
Par_Thruster{2}			= Thrusters{1};		% Ducted main prop
Par_Thruster{2}.r_p		= [-37.8 5 4]';		% thruster position(x,y,z)
Par_Thruster{2}.alpha	= 0;				% thruster orientation

% Thruster 3
% Tunnel thruster, use always positive y-coordinate for defining half breadth of vessel at thruster
Par_Thruster{3}			= Thrusters{5};		% Tunnel thruster
Par_Thruster{3}.r_p		= [35 5 4]';		% thruster position(x,y,z)
Par_Thruster{3}.alpha	= pi/2;				% thruster orientation

% Thruster 4
Par_Thruster{4} 		= Thrusters{5};
Par_Thruster{4}.r_p		= [30 5 4]';	% thruster position(x,y,z)
Par_Thruster{4}.alpha	= pi/2;			% thruster orientation

% Thruster 5
Par_Thruster{5} 		= Thrusters{5};
Par_Thruster{5}.r_p		= [-30 5 4]';	% thruster position(x,y,z)
Par_Thruster{5}.alpha	= pi/2;			% thruster orientation

