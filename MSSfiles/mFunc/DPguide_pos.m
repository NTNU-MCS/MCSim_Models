function [vs,pd,pd_s,vs_s] = DPguide_pos(s,p0,pt,yaw,lambda,lim_vel,lim_acc,p,enable)

% Function that calculates for each s-value, the desired position pd(s),   
% and speed assignment signal vs(t,s), given value of the signal 'enable'. 
% enable = 0: No motion (vs = 0), enable = 1: Motion using sat() function; 
% enable = 2: Motion using smooth tanh() function.
% For path, pd(0) = p0 and pd(1+lambda) = pt.
% The velocity and acceleration are limited according to the 'lim_vel' and 
% 'lim_acc' signals according to the 'p' value. See limR2.m function for 
% details. 
%
% Parameters:
%  s        Non-negative value of path-parameter (typically from integrator).
%  p0       Initial position of path, in R2 [m] (vessel initial position).
%  pt       Target position of path, in R2 [m].
%  yaw      Present vessel yaw angle [radians]
%  lambda   lambda parameter ( note: pd(1+lambda) = pt ).
%  lim_vel  limit velocity of bounding set; scalar for p-norm, R2 vector 
%           for box/diamond constraint.
%  lim_acc  limit acceleration of bounding set; scalar for p-norm, R2 
%           vector for box/diamond constraint.
%  p        norm to use, or 'box' or 'diamond' for box or diamond 
%           constraints.
%  enable   Switch signal that starts (=1) or stops (=0) the motion. 
%  
% Output:
%  vs       Speed assignment signal.
%  pd       Desired path signal.
%  pd_s     Constant partial derivative of path signal with respect to s.
%  vs_s     Partial derivative of speed assignment signal with respect to s.
%
% ____________________________________________________________________
% Author:     Roger Skjetne
% Date:       2019-10-09
% Revisions:  
% ____________________________________________________________________
%

narginchk(9,9);
nargoutchk(2,4);

pd   = ((1+lambda-s)*p0 + s*pt)/(1+lambda);
pd_s = (pt-p0)/(1+lambda);
pd_s_norm = norm(pd_s);

course = atan2(pd_s(2),pd_s(1));
crab   = course - yaw;
vd     = [cos(crab); sin(crab)];

[u_lim,u_vec] = limR2(vd,lim_vel,p,0);
[a_lim,a_vec] = limR2(vd,lim_acc,p,0);

temp = sqrt((1+2*lambda)*(pd_s_norm+eps)*a_lim/2);
if u_lim > temp
    u_lim = temp;
end
kb = (pd_s_norm+eps)/(u_lim^2)*a_lim;

if s < -lambda
    us   = 0.0;
    us_s = 0.0;
elseif s >= -lambda && s < (1-lambda)/2
    if enable > 1.5
        us   = u_lim*tanh(kb*(s+lambda));
        us_s = u_lim*kb*(1-tanh(kb*(s+lambda))^2);
    elseif enable > 0.5
        us   = u_lim*sat(kb*(s+lambda),-1,1);
        if s >= -lambda && s < (1-kb*lambda)/kb
            us_s = u_lim*kb;
        elseif s >= (kb*(1+lambda)-1)/kb && s < (1+kb*(1+lambda))/kb
            us_s = -u_lim*kb;
        else
            us_s = 0.0;
        end
    else
        us   = 0.0;
        us_s = 0.0;
    end
else
    if enable > 1.5
        us   =  u_lim*tanh(kb*(1+lambda-s)); 
        us_s = -u_lim*kb*(1-tanh(kb*(1+lambda-s))^2);
    elseif enable > 0.5
        us   = u_lim*sat(kb*(1+lambda-s),-1,1);
        if s >= -lambda && s < (1-kb*lambda)/kb
            us_s = u_lim*kb;
        elseif s >= (kb*(1+lambda)-1)/kb && s < (1+kb*(1+lambda))/kb
            us_s = -u_lim*kb;
        else
            us_s = 0.0;
        end
    else
        us   = 0.0;
        us_s = 0.0;
    end
end



if enable > 0.5
    vs   = us/(pd_s_norm+eps);
    vs_s = us_s/(pd_s_norm+eps);
else
    vs   = 0.0;
    vs_s = 0;
end

