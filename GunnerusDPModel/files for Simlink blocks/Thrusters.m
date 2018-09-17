% Simplified thruster modelling


% Settings for simplified thruster model, only for Sim.Thrusters == 0

Sim.Thrusters.thrust_sat	= 1;			% 1 for enabling saturation of thrust vector in simple thruster dynamics, 0 for not
Sim.Thrusters.thrust_time	= 1;			% 1 for enabling 1st order thruster dynamics, 0 for not

% Max thrust vector. Only used when Sim.Thrusters.thrust_sat == 1
%Sim.Thrusters.Tau_max		= [10 10 12];

% Time constants for 1st order thruster system, when Sim.Thrusters.thrust_time == 1

Sim.Thrusters.T_surge		= 5;		% s
Sim.Thrusters.T_sway		= 5;		% s
Sim.Thrusters.T_yaw			= 5;		% s
Control.Tau_max             = [1 1 .5] * 10^5; % note, these values should be checked


