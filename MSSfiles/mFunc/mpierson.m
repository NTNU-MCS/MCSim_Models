function [w,S,Tz] = mpierson(Hs,wo,wmax,N)
% [w,S,Tz] = mpierson(Hs,wo,wmax,N) two parameter spectrum for fully developed sea.
% Returns the spectral density function S of the Modified Pierson-Moskowitz 
% spectrum for the frequencies 0 < w < wmax.
%
% Ouputs:
%   w  = vector of equally spaced wave spectrum frequencies (rad/s)
%   S  = vector of power spectral densities (m^2s)
%   Tz = period of zero-crossings (s) 
%
% Inputs:
%   Hs   = significant wave height (m) - mean of the ones third highest waves
%   wo   = peak frequency (rad/s)
%   wmax = maximum wave spectrum frequency:  0 < wo < wmax (rad/s)
%   N    = number of equally spaced frequencies, optionally (default: N=100)
%
% See also pierson, jonswap, torset
%
% Author:   Thor I. Fossen
% Date:     11th August 2001
% Revisions: 

if wo>=wmax, error('wmax must be larger than the peak frequency wo'); end
if nargin==3, N=100; end;
w = 0.01:wmax/N:wmax;

B = 5*wo^4/4;
Tz = (16*pi^3/B)^(1/4);
A = 4*pi^3*Hs^2/Tz^4;

S = A.*w.^(-5).*exp(-B*w.^(-4));

figure(gcf)
plot(w, S/((2*pi/wo)*Hs^2),'g','linewidth',2)
title('Modified Pierson-Moskowitz (PM) spectrum')
xlabel('\omega (rad/s)')
ylabel('Non-dimensional: S(\omega)/(H_s^2T_o)')
