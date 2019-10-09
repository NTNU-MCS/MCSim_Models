function y = sat(x,xmin,xmax)
%
% Returns x saturated between the limits xmin and xmax.
% Parameters:
%  x        Vector, where each element will be saturated.   
%  xmin     Minimum value of x.
%  xmax     Maximum value of x.
%  
% Output:
%  y        Vector of equal dimension as x with saturated values.
%
% ____________________________________________________________________
% Author:     Roger Skjetne
% Date:       2019-04-28
% Revisions:  
% ____________________________________________________________________
%

if xmin > xmax
    disp('Error. xmin is larger than xmax.');
    return;
end

[M,N] = size(x);
y     = zeros(M,N);

for ii=1:M
    for jj=1:N
        if x(ii,jj) < xmin
            y(ii,jj) = xmin;
        elseif x(ii,jj) > xmax
            y(ii,jj) = xmax;
        else
            y(ii,jj) = x(ii,jj);
        end
    end
end

