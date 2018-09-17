function angle = rad2pipi(angle);
% RAD2PIPI Converts an angle in rad to the interval <-pi pi>
%          Should be applied to all heading errors in a feedback control system
%          in order to avoid discontinuities.
%
% Author:   Thor I. Fossen
% Date:     21th July June 2001
% Revisions: 

% convert angle in rad to the interval <-pi pi> 
angle = angle - fix(angle/(2*pi))*2*pi;
angle = angle - fix(angle/pi)*2*pi;