% ###########################################################################################################
% 
% File: "WF_bodep.m"
%
% Project: "Extreme sea observer and controller design for MCSim"
%    
% Abstract:     Displays the Bode plot of the transfer function of the nonlinear passive observer.
%               Demonstrates the notch effect from wave filtering.
%               
% Inputs:       Control.VessObs.lmb         - Assumed wave spectrum damping from "dysystem.m"
%               Control.VessObs.omg_0       - Calculated wave peak-frequency from MCSim
%               Control.VessObs.K1          - Calculated filter gain K1 from "dpsystem.m"
%               Control.VessObs.K2          - Calculated filter gain K2 from "dpsystem.m"
%               Control.VessObs.K3          - Calculated filter gain K3 from "dpsystem.m"
%               Control.VessObs.K4          - Calculated filter gain K4 from "dpsystem.m"
%
% Outputs:      -
%
% Calls:        -
%
% Author:       Martin Austad, Autumn 2003
%
% Modified:     -
%
% ########################################################################################################### 

% Renaming variables, for the sake of simplified overview in calculations
%Control.VessObs.damp  = 1;
%Control.VessObs.omg_c = 1.22*Control.VessObs.omg_01;
%Control.VessObs.k     = 1.22;

lambda = .1;%Control.VessObs.lmb;
wo     = Control.VessObs.omg_0{1};

%Control.VessObs.K41 = [.1 0 0; 0 .18 0; 0 0 .72];  % Dong
%Control.VessObs.K31 = diag([.01 1000 50000]);

%Control.VessObs.K41 = [20000 0 0; 0 .18 0; 0 0 .72];  % Dong
%Control.VessObs.K31 = diag([200 1000 50000]);

%K1     = Control.VessObs.K11;
%K2     = Control.VessObs.K21;
%K3     = Control.VessObs.K31;
%K4     = Control.VessObs.K41;
K1     = [zeros(3);zeros(3)];
K2     = Control.VessObs.K2{1};
K3     = Control.VessObs.K3{1};
K4     = Control.VessObs.K4{1};


i=1;  % SURGE MOTION

h0 = tf([1 2*lambda*wo wo^2], [1 (K1(i+3,i) + K2(i,i) + 2*lambda*wo) (wo^2+2*lambda*wo*K2(i,i) -K1(i,i)*wo^2) (wo^2*K2(i,i))]);
     
hB = tf(K4(i,i)*[1 (K3(i,i)/K4(i,i))],[1 Control.VessObs.T_inv(i,i)]);        
wave = tf([wo^2 0],[1 2*lambda*wo  wo^2]);

% BODE PLOT 
% -----------------------------------------------------------------------------------------------------------
figure

w = logspace(-4,1.5,100);
bode(series(h0,hB),wave,w);
grid on
% -----------------------------------------------------------------------------------------------------------


i=2;  % SWAY MOTION

h0 = tf([1 2*lambda*wo wo^2], [1 (K1(i+3,i) + K2(i,i) + 2*lambda*wo) (wo^2+2*lambda*wo*K2(i,i) -K1(i,i)*wo^2) (wo^2*K2(i,i))]);
    
hB = tf(K4(i,i)*[1 (K3(i,i)/K4(i,i))],[1 Control.VessObs.T_inv(i,i)]);        
wave = tf([wo^2 0],[1 2*lambda*wo  wo^2]);
 
% BODE PLOT 
% -----------------------------------------------------------------------------------------------------------
figure
figure(gcf)
bode(series(h0,hB),wave,w);
grid on
% -----------------------------------------------------------------------------------------------------------


i=3;  % YAW MOTION

h0 = tf([1 2*lambda*wo wo^2], [1 (K1(i+3,i) + K2(i,i) + 2*lambda*wo) (wo^2+2*lambda*wo*K2(i,i) -K1(i,i)*wo^2) (wo^2*K2(i,i))]);
     
hB = tf(K4(i,i)*[1 (K3(i,i)/K4(i,i))],[1 Control.VessObs.T_inv(i,i)]);        
wave = tf([wo^2 0],[1 2*lambda*wo  wo^2]);
 
% BODE PLOT 
% -----------------------------------------------------------------------------------------------------------
figure
figure(gcf)

bode(series(h0,hB),wave,w);
grid on
% -----------------------------------------------------------------------------------------------------------

% ###########################################################################################################
% END OF FILE.
% ###########################################################################################################