%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
%
% MCS.m
%
% Project:	"Marine Cybernetics Simulator: MCSim"
%
% Abstract:	Run at startup. Reset screen and set paths. Display welcome message
% 
% Inputs:	-
%
% Outputs:	-
%	
% Calls:	-
%		
% Author:	Øyvind Smogeli, Summer 2003
%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 

% Clear window
clc;


% Save old path to file
oldpath = path;
save Initialize\oldpath_file oldpath

% Set paths for the subfolders, independent of working directory
pathstring = which('MCS');
pathstring = pathstring(1:(length(pathstring)-6));

path(path, pathstring);
path(path, strcat(pathstring,'\Process'));
path(path, strcat(pathstring,'\Initialize'));
path(path, strcat(pathstring,'\Vessel_data'));
path(path, strcat(pathstring,'\Results'));
path(path, strcat(pathstring,'\Toolbox'));


% Display welcome message
disp(' ')
disp('***********************************************************')
disp('*              WELCOME TO MCSIM Beta 4                    *')
disp('***********************************************************')
disp(' ')
disp('   - Set simulation parameters in "Parameters.m" ')
disp('   - Set controller parameters in "Controls.m" ')
disp('   - Define thruster configuration in "Thruster_config.m" ')
disp('   - Initialize simulation with "Init.m" ')
disp('   - Run simulator from Simulink file "MCSim.mdl" ')
disp('   - Save results with "Savedata.m" ')
disp('   - Plot results with "Vessel_Results.m", "Control_Results.m"')
disp('     and "Thruster_Results.m" ')
disp('   - Check additional info in "Readme.m" ')
disp(' ')
disp('   - GNC toolbox must be installed, get it at:')
disp('     http://www.marinecybernetics.com/software.htm ')
disp(' ')
disp('   - Subfolder paths set. To reset old path, run script "reset_path.m".')
disp(' ')
disp('***********************************************************')

