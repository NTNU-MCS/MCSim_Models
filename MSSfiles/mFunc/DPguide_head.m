function [vs,psid,psid_s,vs_s] = DPguide_head(s,psi0,psit,lambda,lim_yawrate,lim_yawacc,enable)

% Function that calculates for each s-value, the desired heading psid(s),   
% and speed assignment signal vs(t,s), given value of the signal 'enable'. 
% enable = 0: No yaw motion (vs = 0); enable = 1: Motion using sat() 
% function;. enable = 2: Motion using smooth tanh() function.
% For heading, psid(0) = psi0 and psid(1+lambda) = psit.
% The yaw rate and acceleration are limited according to the 'lim_yawrate' 
% and 'lim_yawacc' parameters.  
%
% Parameters:
%  s        Non-negative value of path-parameter (typically from integrator).
%  psi0     Initial position of path, in R2 [m] (vessel initial position).
%  psit     Target position of path, in R2 [m].
%  lambda   lambda parameter ( note: pd(1+lambda) = pt ).
%  lim_yawrate  limit yaw rate value.
%  lim_yawacc  limit yaw acceleration value.
%  enable   Switch signal that starts (=1) or stops (=0) the motion. 
%  
% Output:
%  vs       Speed assignment signal.
%  psid     Desired path signal.
%  psid_s   Constant partial derivative of path signal with respect to s.
%  vs_s     Partial derivative of speed assignment signal with respect to s.
%
% ____________________________________________________________________
% Author:     Roger Skjetne
% Date:       2019-10-09
% Revisions:  
% ____________________________________________________________________
%

narginchk(7,7);
nargoutchk(2,4);

psid   = ((1+lambda-s)*psi0 + s*psit)/(1+lambda);
psid_s = (psit-psi0)/(1+lambda);
psid_s_norm = norm(psid_s);

temp = sqrt((1+2*lambda)*(psid_s_norm+eps)*lim_yawacc/2);
if lim_yawrate > temp
    lim_yawrate = temp;
end
kr = (psid_s_norm+eps)/(lim_yawrate^2)*lim_yawacc;

if s < -lambda
    rs   = 0.0;
    rs_s = 0.0;
elseif s >= -lambda && s < (1-lambda)/2
    if enable > 1.5
        rs   = lim_yawrate*tanh(kr*(s+lambda));
        rs_s = lim_yawrate*kr*(1-tanh(kr*(s+lambda))^2);
    elseif enable > 0.5
        rs   = lim_yawrate*sat(kr*(s+lambda),-1,1);
        if s >= -lambda && s < (1-kr*lambda)/kr
            rs_s = lim_yawrate*kr;
        elseif s >= (kr*(1+lambda)-1)/kr && s < (1+kr*(1+lambda))/kr
            rs_s = -lim_yawrate*kr;
        else
            rs_s = 0.0;
        end
    else
        rs   = 0.0;
        rs_s = 0.0;
    end
else
    if enable > 1.5
        rs   =  lim_yawrate*tanh(kr*(1+lambda-s)); 
        rs_s = -lim_yawrate*kr*(1-tanh(kr*(1+lambda-s))^2);
    elseif enable > 0.5
        rs   = lim_yawrate*sat(kr*(1+lambda-s),-1,1);
        if s >= -lambda && s < (1-kr*lambda)/kr
            rs_s = lim_yawrate*kr;
        elseif s >= (kr*(1+lambda)-1)/kr && s < (1+kr*(1+lambda))/kr
            rs_s = -lim_yawrate*kr;
        else
            rs_s = 0.0;
        end
    else
        rs   = 0.0;
        rs_s = 0.0;
    end
end



if enable > 0.5
    vs   = rs/(psid_s_norm+eps);
    vs_s = rs_s/(psid_s_norm+eps);
else
    vs   = 0.0;
    vs_s = 0;
end

