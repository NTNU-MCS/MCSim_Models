% Hydrodynamic processing of ShipX RV Gunnerus data
% vessel = veres2vessel('input.out');
% vesselABC = vessel2ss(vessel);  % fluid memory effects
% 
% plot
% plotABC(vessel,'A')
% plotABC(vessel,'B')
% plotTF(vessel,'motion','rads',1)
% plotTF(vessel,'force','rads',1)
% plotWD(vessel,'rads',1)
% 
% display main data
% display(vessel.main);


vessel = veres2vessel('GunnerusMan_input');  % reads .re1, .re2, .re7, .re8 and .hyd
vesselABC = vessel2ssGunnerus(vessel);  % calculated fluid memory effects, A(0), B(0)




% read_veres_TF(strcat('input','.re1'), 0);
% Outputs:
%
% TF_data contains transfer function data extracted from the file in form
% of amplitude amplification and phase, structured as:
%
% 	vessel.forceRAO.amp{dofno}(freqno,headnno,velno)   
% 	vessel.forceRAO.phase{dofno}(freqno,headnno,velno) 
%
% 	vessel.motionRAO.amp{dofno}(freqno,headnno,velno)  
% 	vessel.motionRAO.phase{dofno}(freqno,headnno,velno) 
% 	vessel.motionRAO.w(1,freqno)
%
% 	vessel.headings   	
% 	vessel.velocities	
% 
% 'vel', 'dir' and 'freq' contain velocities in m/s, directions in radians
% and frequencies in rad/s respectively.
%
% 'Amp' and 'Phase' are 3D matrices with dimensions (freq, phase, vel)

