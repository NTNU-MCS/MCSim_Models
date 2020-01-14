function [y,v] = limR2(x,lim,p,plotflag)

% Returns the limiting value in the direction of an R2 vector according 
% to a specified p-norm constraing, 'box' constraint, or 'diamond'
% constraint. For p-norms it is assumed symmetry in all directions. For box
% or diamond constraints, it is assumed port-starboard and fore-aft
% symmetries.
%
% Parameters:
%  x        2xN matrix, where each column is a vector in R2 to be bounded.   
%  lim      limit value of bounding set; scalar for p-norm, R2 vector for 
%           box/diamond constraint.
%  p        norm to use, or 'box' or 'diamond' for box or diamond 
%           constraints.
%  plotflag 1: plotting x and y; 0: No plot. No plot if omitted.
%  
% Output:
%  y        Vector of equal dimension as each row of x with limit values.
%  v        2xN matrix, where each column is the saturated vector.
%
% ____________________________________________________________________
% Author:     Roger Skjetne
% Date:       2019-05-02
% Revisions:  
% ____________________________________________________________________
%

narginchk(3,4);
nargoutchk(1,2);

type = -1;

if ~isnumeric(x)
    error('x must be a vector or matrix of numbers. Should be R2, that is, two rows if a matrix of vectors.');
else
    [M,N] = size(x);
    if M ~= 2
        error('x has wrong dimension. Should be R2, that is, two rows if a matrix of vectors.');
    end
    beta = atan2(x(2,:),x(1,:));
end
if ~isnumeric(lim)
    error('lim must be a scalar or R2 vector of positive numbers. ');
else
    if length(lim)>2 
        error('lim must be a scalar or R2 vector of positive numbers. ');
    elseif length(lim)==1 && lim<=0
        error('lim must be a scalar or R2 vector of positive numbers. ');
    elseif length(lim)==2 && (lim(1)<=0 || lim(2)<=0)
        error('lim must be a scalar or R2 vector of positive numbers. ');
    end
end
if isnumeric(p)
    if p<1
        error('p must be a number (p-norm) within [1,inf], or it must be a string with "box" or "diamond".');
    else
        type = 1; % p correspond to p-norm constraint
    end
elseif ischar(p)
    if strcmp(p,'box')
        type = 2; % p correspond to box constraint
        if length(lim)~=2
            error('lim must be an R2 vector of numbers if p = "box" or "diamond. ');
        end
    elseif strcmp(p,'diamond')
        type = 3; % p correspond to diamond constraint
        if length(lim)~=2
            error('lim must be an R2 vector of numbers if p = "box" or "diamond". ');
        end
    else
        error('p must be a number (p-norm) within [1,inf], or it must be a string with "box" or "diamond".');
    end
else
    error('p must be a number (p-norm) within [1,inf], or it must be a string with "box" or "diamond".');
end

y = zeros(1,N);
v = zeros(M,N);

switch type
    case 1  % p-norm
        for k=1:N
            temp = lim*x(:,k)/norm(x(:,k),p);
            y(1,k) = norm(temp);
            if y(1,k) <= norm(x(:,k))
                v(:,k) = temp;
            else
                v(:,k) = x(:,k);
            end

        end
        
    case 2  % box
        for k=1:N
%            beta = angle(x(1)+x(2)i);
           if cos(beta(k))~= 0
                v1 = [lim(1); lim(1)*tan(beta(k))];
                v2 = [lim(2)/tan(beta(k)); lim(2)];
                magn = min(norm(v1),norm(v2));
                y(1,k) = magn;
                if magn <= norm(x(:,k))
                    v(:,k) = magn*[cos(beta(k));sin(beta(k))];
                else
                    v(:,k) = x(:,k);
                end
%            elseif cos(beta)<0

           else
                magn = abs(lim(2));
                v(:,k) = magn*[cos(beta(k));sin(beta(k))];
                y(1,k) = magn;
           end
        end
        
    case 3  % diamond
        for k=1:N
            if x(1,k) >= 0 && x(2,k) >= 0       % 1st quadrant
                p0 = [0;lim(2)]; 
                p1 = [lim(1);0];
            elseif x(1,k) < 0 && x(2,k) >= 0    % 2nd quadrant
                p0 = [-lim(1);0]; 
                p1 = [0;lim(2)];
            elseif x(1,k) < 0 && x(2,k) < 0     % 3rd quadrant
                p0 = [-lim(1);0]; 
                p1 = [0;-lim(2)];
            else                                % 4th quadrant
                p0 = [0;-lim(2)]; 
                p1 = [lim(1);0];
            end
            ts   = [x(:,k) p0-p1]\p0;
            temp = ts(1)*x(:,k);
            magn = norm(temp);
            if magn <= norm(x(:,k))
                v(:,k) = temp;
            else
                v(:,k) = x(:,k);
            end
            y(1,k) = magn;
        end
        
    otherwise
end

    
if plotflag
    x_lim = y.*[cos(beta);sin(beta)];
    plot(x(2,:),x(1,:),'r',x_lim(2,:),x_lim(1,:),'b',v(2,:),v(1,:),'g');
    legend('x','limit','v'); grid; axis equal;
end