%% Dynamical system configuration File
% Project description 
System.Project          = 'Example: Dynamical System Simulation';
System.Model            = 'Test model';


%% System parameters
System.SampleTime       = -1;
System.NumContStates	= 3;
System.NumDiscStates	= 0;
System.NumInputSignals	= 1;
System.NumOutputSignals	= 4;
System.NumParameters	= 3;
System.ConfigFile       = 'sDynSys_param';


%% Matlab calculations:
k1 = 2; k2 = 1;
A = [0 1; -k1 -k2];
n = 2;
Q = [10*(n-1)+1.5 0.7; 0.7 .6];
P = lyap(A',Q);


%% Model parameters
System.Param{1}.tag  	= 'A';
System.Param{1}.values	= A;

System.Param{2}.tag  	= 'mu';
System.Param{2}.values	= .03;

System.Param{3}.tag  	= 'P';
System.Param{3}.size	= 4;
System.Param{3}.values	= P;
% System.Param{3}.values	= 0.5*eye(2);


[status, error] = m2ini(System, System.ConfigFile);

%% Model mask parameters
SimParam = [System.SampleTime, System.NumContStates, System.NumDiscStates, ...
            System.NumInputSignals, System.NumOutputSignals, System.NumParameters]';
Fname    = [System.ConfigFile,'.ini'];
Xc0      = [-sqrt(2)/2 sqrt(2)/2 pi]';
Xd0      = [];