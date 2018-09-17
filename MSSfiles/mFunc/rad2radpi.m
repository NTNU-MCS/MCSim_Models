% Program that maps an angle in radians into [-pi, pi).
% Roger Skjetne - 05.09.2003
function y = rad2radpi(x);

r = rem(x+sign(x)*pi,2*pi);
y = r - sign(x)*pi;