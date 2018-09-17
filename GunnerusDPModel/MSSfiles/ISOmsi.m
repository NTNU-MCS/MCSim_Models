function [a_z,w_e] = ISOmsi(t)
% msi = ISOMSI(a_z,w_e) computes the Motion Sickness Incidence using the method 
%       of O'Hanlon and McCauley (1974).
%
% a_z: mean of absolute vertical acceleration (m/s^2), i.e. a_z = mean(abs(a_measured(:,1)))
% w_e: enconter frequency (rad/s), see encounter.m.
% msi: the percentage of persons that become seasick during a 2 hours sail
%
% Refs.  - A. R. J. M. Lloyd (1989). Seakeeping Behaviour in Rough Water. Ellis Horwoowd Ltd.
%        - E. V. Lewis (Ed.) (1989). Principles of Naval Architecture. Vol III Motions in
%             Waves and Controllability, 2nd. ed., SNAME. 
%        - J.F. O'Hanlon and M. E. McCauley (1974). Motion Sickness Incidence as a Function of
%             Vertical Sinusoidal Motion. Aerospace Medicine AM-45(4):366-369.
%
% Author:    Thor I. Fossen
% Date:      5th November 2001
% Revisions: 

% ISO 2631-3, 1997 Motion Sickness Index (MSI)

i = 1;
for f = 0.1:0.03:0.63,
    w_e(i) = f*2*pi;
    if (0.1 <= f & f < 0.315),
        a_z(i) = 0.5*sqrt(2/t);
    else
        a_z(i) = 0.5*sqrt(2/t)*6.8837*f^1.67;
    end
    i = i+1;
end