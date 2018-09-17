function [w,S] = torset(Hs,wo,wmax,N)
% [w,S] = torset(Hs,wo,wmax,N) is an empirical two peaked spectrum for swell and 
% developing sea based on experimental data from the North Sea. For small peak
% frequencies, i.e.  0 < wmax <= 0.6 only one peak in the spectrum appears.
% Returns the spectral density function S of the Torsethaugen spectrum for 
% the frequencies:  0 < w < wmax (rad/s).
%
% Ouputs:
%   w     = vector of equally spaced wave spectrum frequencies (rad/s)
%   S     = vector of power spectral densities (m^2s)
%
% Inputs:
%   Hs    = significant wave height (m) - mean of the ones third highest waves
%   wo    = peak frequency (rad/s)
%   wmax  = maximum wave spectrum frequency:  0 < wo < wmax (rad/s)
%   N     = number of equally spaced frequencies, optionally (default: N=100)
%
% See also pierson, mpierson, jonswap
%
% Ref: K.Torsethaugen (1996): "Model for a Doubly Peaked Wave Spectra"
%      Sintef report no.: STF22 A96204 prepared for Norsk Hydro.
%
% Author:    G. Kleiven, Norsk Hydro 
% Date:      15th May 2000
% Revisions:  6th June 2001, Svein I. Sagatun, Norsk Hydro - minor revisions
%            14th August 2001, Thor I. Fossen - IO compatibel with the GNC toolbox

%-----------------------------------------------------------------------------
% IO for Matlab GNC toolbox, T. I. Fossen
%-----------------------------------------------------------------------------
if wo>=wmax, error('wmax must be larger than the peak frequency wo'); end
if wmax<=0.6, disp('choose wmax>0.60 to obtain a two peaked spectrum'); end
if nargin==3, N=100; end

Hmo  = Hs;           % significant wave height (m)
Tp   = 2*pi/wo;      % peak period (s)
fmax = wmax/(2*pi);  % maximum frequency (Hz)
Nfrq = N;            % number of frequencies
%---------------------------------------------------------------------------
% source code: Norsk Hydro
%---------------------------------------------------------------------------
f2pii = 2*pi;
domg = f2pii*fmax/Nfrq;
fwtp = f2pii/Tp;
omg = [1:Nfrq]*domg;
%---------------------------------------------------------------------------
%  Parameters:
%---------------------------------------------------------------------------
af = 6.6;
ae = 2.0; 
au = 25;   
a10 = 0.7;
a1 = 0.5;
kg = 35.0;
kg0 = 3.5;
kg1 = 1.0;
r = 0.857;
k0 = 0.5;
k00 = 3.2;
m0 = 4.0;
b1 = 2.0;
a20 = 0.6;
a2 = 0.3;
a3 = 6.0;
s0 = 0.08;
s1 = 3.0;
b2 = 0.7;
b3 = 3.0;
sigma_a2 = 2*0.07^2;
sigma_b2 = 2*0.09^2;     
tf = af*Hmo^(1./3);
      
if Tp < tf;     
%-------------------------------------------------------------------------
% Predominant wind sea peak:
%-------------------------------------------------------------------------
	 tl = ae*(Hmo^0.5);
	 eps1 = (tf-Tp)/(tf-tl);
	 rpw = (1-a10)*exp(-((eps1/a1)^2))+a10;
	 hsw = rpw*Hmo;
	 hss = sqrt(1.-rpw^2)*Hmo;
	 tpw = Tp;
	 tps = tf+b1;
	 sp = ((f2pii/9.81)*hsw/(Tp^2));
	 gammaw = kg*(1+kg0*exp(-Hmo/kg1))*(sp^r);
	 gammas = 1.;
	 nw = k0*sqrt(Hmo)+k00;
	 mw = m0;
	 ns = nw;
	 ms = mw;
    if ms < 1
       ms = 1;
    end         % Constraint implemented by GKl
	 g_argw = (nw-1)/mw;
	 g_args = (ns-1)/ms;
	 if g_args < 0
	    g0s = 1./((1./ms)*gamma(g_args)*((ns/ms)^(-g_args)));
	 else
	    g0s = 1./((1./ms)*gamma(g_args)/((ns/ms)^(g_args)));
	 end
	 if g_argw < 0
	    g0w = 1./((1./mw)*gamma(g_argw)*((nw/mw)^(-g_argw)));
	 else
	    g0w = 1./((1./mw)*gamma(g_argw)/((nw/mw)^(g_argw)));
	 end
	 a1m = 4.1;
	 b1m = 2.0*(mw^0.28)-5.3;
	 c1m = -1.45*(mw^0.1)+0.96;
	 a2m = 2.2/(mw^3.3)+0.57;
	 b2m = -0.58*mw^0.37+0.53;
	 c2m = -1.04/(mw^1.9)+0.94;
	 if c1m < 0 
	    f1w = a1m / (nw-b1m)^(-c1m);
	 else
	    f1w = a1m * (nw-b1m)^c1m;
	 end
	 if b2m < 0
	    f2w = a2m / nw^(-b2m) + c2m;
	 else
	    f2w = a2m * nw^(b2m) + c2m;
	 end
	 b1m = 2.0*(ms^0.28)-5.3;
	 c1m = -1.45*(ms^0.1)+0.96;
	 a2m = 2.2/(ms^3.3)+0.57;
	 b2m = -0.58*ms^0.37+0.53;
	 c2m = -1.04/(ms^1.9)+0.94;
	 if c1m < 0 
	    f1s = a1m / (ns-b1m)^(-c1m);
	 else
	    f1s = a1m * (ns-b1m)^c1m;
	 end
	 if b2m < 0
	    f2s = a2m / ns^(-b2m) + c2m;
	 else
	    f2s = a2m * ns^(b2m) + c2m;
	 end
	 
	 agammaw = (1+f1w*log(gammaw)^(f2w))/gammaw;
	 agammas = (1+f1s*log(gammas)^(f2s))/gammas;
else
%----------------------------------------------------------------------------------------------------------------
% Predominant swell peak:
%--------------------------------------------------------------------------
	 tu = au;
	 epsu = (Tp-tf)/(tu-tf);
	 rps = (1.-a20)*exp(-(epsu/a2)^2)+a20;
	 hss = rps*Hmo;
	 hsw = sqrt(1.-rps^2)*Hmo;
	 tps = Tp;
	 ns = k0*sqrt(Hmo)+k00;
	 ms = m0;
	 nw = ns;
	 mw = m0*(1-b2*exp(-Hmo/b3));
	 s4 = s0*(1-exp(-Hmo/s1));
	 g_argw = (nw-1)/mw;
	 g_args = (ns-1)/ms;
	 if g_args < 0
       g0s = 1./((1./ms)*gamma(g_args)*((ns/ms)^(-g_args)));
    else
	    g0s = 1./((1./ms)*gamma(g_args)/((ns/ms)^(g_args)));
	 end
	 if g_argw < 0
	    g0w = 1./((1./mw)*gamma(g_argw)*((nw/mw)^(-g_argw)));
	 else
	    g0w = 1./((1./mw)*gamma(g_argw)/((nw/mw)^(g_argw)));
	 end
	 tpw = ((g0w*hsw^2)/(16*s4*(0.4^nw)))^(1./(nw-1.));
	 sf = ((f2pii/9.81)*Hmo/(tf^2));
	 gammaw = 1.;
	 gamma_f = kg*(1+kg0*exp(-Hmo/kg1))*sf^r;
	 gammas = gamma_f*(1.+a3*epsu);
	 a1m = 4.1;
	 b1m = 2.0*(mw^0.28)-5.3;
	 c1m = -1.45*(mw^0.1)+0.96;
	 a2m = 2.2/(mw^3.3)+0.57;
	 b2m = -0.58*(mw^0.37)+0.53;
	 c2m = -1.04/(mw^1.9)+0.94;
	 if c1m < 0
	    f1w = a1m / (nw-b1m)^(-c1m);
	 else
	    f1w = a1m * (nw-b1m)^c1m;
	 end
	 if b2m < 0
	    f2w = a2m / nw^(-b2m) + c2m;
	 else
	    f2w = a2m * nw^(b2m) + c2m;
	 end
	 b1m = 2.0*(ms^0.28)-5.3;
	 c1m = -1.45*(ms^0.1)+0.96;
	 a2m = 2.2/(ms^3.3)+0.57;
	 b2m = -0.58*(ms^0.37)+0.53;
	 c2m = -1.04/(ms^1.9)+0.94;
	 if c1m < 0
	    f1s = a1m / (ns-b1m)^(-c1m);
	 else
	    f1s = a1m * (ns-b1m)^c1m;
	 end
	 if b2m < 0
	    f2s = a2m / ns^(-b2m) + c2m;
	 else
	    f2s = a2m * ns^(b2m) + c2m;
	 end
	 agammaw = (1+f1w*log(gammaw)^(f2w))/gammaw;
	 agammas = (1+f1s*log(gammas)^(f2s))/gammas;
end
      
fdenorm_s = (tps*(hss^2))/16;
fdenorm_w = (tpw*(hsw^2))/16;

%===================================================================
%  Estimates spectral density for each frequency in array omg:
%===================================================================

f = omg/f2pii;
%-------------------------------------------------------------------
%Wind sea contribution:
%-------------------------------------------------------------------
fnw = f*tpw;
in = max(find(fnw < 1));
ftest1(1:in) = exp(-(((fnw(1:in)-1).^2)/sigma_a2));
ftest1(in+1:Nfrq) = exp(-(((fnw(in+1:Nfrq)-1).^2)/sigma_b2));
gamma_wf = gammaw.^ftest1;
gamma_ws_1 = fnw.^(-nw);
gamma_ws_2 = exp(-(nw/mw)*fnw.^(-mw));
gamma_ws = gamma_ws_1.*gamma_ws_2;
sw = g0w*agammaw*gamma_ws.*gamma_wf*fdenorm_w;
%----------------------------------------------------------------------
%Swell contribution:
%----------------------------------------------------------------------
fns = f*tps;
is = max(find(fns < 1));
ftest2(1:is) = exp(-(((fns(1:is)-1).^2)/sigma_a2));
ftest2(is+1:Nfrq) = exp(-(((fns(is+1:Nfrq)-1).^2)/sigma_b2));
gamma_sf = gammas.^ftest2;
gamma_ss_1 = fns.^(-ns);
gamma_ss_2 = exp(-(ns/ms)*fns.^(-ms));
gamma_ss = gamma_ss_1.*gamma_ss_2;
ss = g0s*agammas*gamma_ss.*gamma_sf*fdenorm_s;
%-----------------------------------------------------------------------------
%Estimates spectral-density Sf(m^2*s)
%-----------------------------------------------------------------------------
sf = sw + ss;  % frq. in Hz

l = max(size(sf));
df = ones(l,1)*domg/f2pii;
var_sf = sf*df;
Hmo_est = 4*sqrt(var_sf);

%-----------------------------------------------------------------------------
% IO for Matlab GNC toolbox, T. I. Fossen
%-----------------------------------------------------------------------------
w = 2*pi*f;
S = sf;

figure(gcf)
plot(w, S/((2*pi/wo)*Hs^2),'b','linewidth',2)
title('Torsethaugen spectrum')
xlabel('\omega (rad/s)')
ylabel('Non-dimensional: S(\omega)/(H_s^2T_o)')