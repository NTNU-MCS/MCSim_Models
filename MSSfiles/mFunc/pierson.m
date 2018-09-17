function [w,S,wo,V] = pierson(Hs,wmax,N)
% [w,S,wo,V] = pierson(Hs,wmax,N) computes the Pierson-Moskowitz (PM) 
% power spectral density function S for the frequencies 0 < w < wmax.
%
% Ouputs:
%   w  = vector of equally spaced wave spectrum frequencies (rad/s)
%   S  = vector of power spectral densities (m^2s)
%   wo = peak frequency (rad/s)
%   V  = wind speed (m/s)
%
% Inputs:
%   Hs   = significant wave height (m) - mean of the ones third highest waves
%   wmax = maximum wave spectrum frequency:  0 < wo < wmax (rad/s)
%   N    = number of equally spaced frequencies, optionally (default: N=100)
%
% See also mpierson, jonswap, torset
%
% Author:   Thor I. Fossen
% Date:     11th August 2001
% Revisions: 

if nargin==2, N=100; end
w = 0.01:wmax/N:wmax;

g = 9.81;
A = 0.0081*g^2;
B = 3.11/(Hs^2);
S = A.*w.^(-5).*exp(-B*w.^(-4));
wo = (4*B/5)^(1/4);
V = sqrt((Hs*g/0.21));

figure(gcf)
plot(w, S/((2*pi/wo)*Hs^2),'c','linewidth',2)
title('Pierson-Moskowitz spectrum')
xlabel('\omega (rad/s)')
ylabel('Non-dimensional: S(\omega)/(H_s^2T_o)')
