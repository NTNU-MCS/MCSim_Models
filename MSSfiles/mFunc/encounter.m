function w_e = encounter(w_o,U,beta)
% w_e = ENCOUNTER(w_o,U,beta) computes the encounter frequency w_e (rad/s) as a function
%       of wave peak frequeny w_o (rad/s), vessel speed U (m/s) and wave direction 
%       beta, 0 (deg) for following seas and 180 (deg) for head seas.
%
% Author:   Thor I. Fossen
% Date:     2nd November 2001
% Revisions: 

g = 9.8; % acceleration of gravity
w_e = abs(w_o - w_o.^2*U*cos(beta*pi/180)/g);