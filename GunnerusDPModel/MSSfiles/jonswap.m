function [w,S,Tz] = jonswap(Hs,wo,wmax,gamma,N)
% [w,S,Tz] = JONSWAP(Hs,wo,wmax,gamma,N) two parameter spectrum for developing sea.
% Returns the spectral density function S of the JONSWAP spectrum for the 
% frequencies:  0 < w < wmax (rad/s).
%
% Ouputs:
%   w  = vector of equally spaced wave spectrum frequencies (rad/s)
%   S  = vector of power spectral densities (m^2s)
%   Tz = period of zero-crossings (s) 
%
% Inputs:
%   Hs    = significant wave height (m) - mean of the ones third highest waves
%   wo    = peak frequency (rad/s)
%   wmax = maximum wave spectrum frequency:  0 < wo < wmax (rad/s)
%   gamma = spectrum peak factor (optionally, default gamma = 3.3)
%   N     = number of equally spaced frequencies, optionally (default: N=100)
%
% See also pierson, mpierson, torset
%
% Author:   Thor I. Fossen
% Date:     14th August 2001
% Revisions: 

if wo>=wmax, error('wmax must be larger than the peak frequency wo'); end
if nargin==3, gamma = 3.3; N= 100; end  
if nargin==4, N = 100; end  
w = 0.01:wmax/N:wmax;

% wo = (4*B/5)^(1/4) and B = 944/T1^4;
B = 5*wo^4/4;
T1 = (944/B)^(1/4);
A = 155*(Hs^2/T1^4);

for i=1:N,
    if w(i)>5.24/T1,
        sigma(i) = 0.09;
    else
        sigma(i) = 0.07;
    end
end

Y = exp(-((0.191.*w*T1-1)./(sqrt(2).*sigma)).^2);
S = A.*w.^(-5).*exp(-B*w.^(-4)).*gamma.^Y;
Tz = T1/1.073;

figure(gcf)
plot(w, S/((2*pi/wo)*Hs^2),'r','linewidth',2)
title('JONSWAP spectrum')
xlabel('\omega (rad/s)')
ylabel('Non-dimensional: S(\omega)/(H_s^2T_o)')

