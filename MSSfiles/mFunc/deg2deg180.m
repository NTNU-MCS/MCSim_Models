% Program that maps an angle in degrees into [-180, 180).
% Roger Skjetne - 05.09.2003
function y = deg2deg180(x);

r = rem(x+sign(x)*180,360);
y = r - sign(x)*180;