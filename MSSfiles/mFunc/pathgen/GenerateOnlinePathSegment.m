    function Path = GenerateOnlinePathSegment(p0,pt,t0,r,lambda,PlotHandle)
%
% Path = GenerateOnlinePathSegment(p0,pt,t0,r,lambda,PlotFlag)
%
% Function that generates the coefficients for a path segment from initial p0 to target pt.
% Each subpath is parametrized by t in [0,1), and smoothly connected at the waypoints. This 
% then correspond to a hybrid parametrized path, where Index i identifies the subpath, and 
% the parameter t in [0,1) identifies the location along the subpath. The path and its 
% derivatives up to r = (Order-1)/2 are guaranteed to be continuous at the waypoints.
%
% Input data:
%
%    p0         - A vector [x0,y0] of initial waypoint.
%    pt         - A vector [xt,yt] of target waypoint.
%    t0         - A vector [tx0,ty0] of tangent vector of path at initial waypoint p0. 
%    r          - Order of derivatives with continuity at the waypoints. 
%    lambda     - Constant that tunes the curvature at p0 and pt. 
%    PlotHandle - Figure number for plotting the path; Empty: No plot.
%
% Output data:
%
%    Path       - A Matlab data structure with the following fields:
%                .NumSubpaths: Number of subpaths.
%                .Order: Order of polynomials (=2*ord+1).
%                .WP: x- and y-coordinates of waypoints.
%                .LinSys: Liner set of equations Ax=b to solve the subpaths:
%                   .A:     Common A matrix for both x- and y-coefficients.
%                   .bx:    Cell structure that contains the b-vector for each subpath for the x-coordinates.
%                   .by:    Cell structure that contains the b-vector for each subpath for the y-coordinates.
%                .coeff: Coefficients for the subpaths:
%                   .a:     Cell structure that contains the a-coefficients for the x-coordinates.
%                   .b:     Cell structure that contains the b-coefficients for the y-coordinates.
%                   .a_der: Cell structure that contains the a-coefficients for the derivatives 
%                           for the x-coordinates of the path up to the polynomial order.
%                   .b_der: Cell structure that contains the b-coefficients for the derivatives 
%                           for the y-coordinates of the path up to the polynomial order.
%
%
%    Copyright: 	Roger Skjetne, NTNU
%    Author:        Roger Skjetne
%    Date created:  2019.03.02  Roger Skjetne.
%    Revised:      	
%


%% Initialization
ord  = 2*r+1;
% if mod(ord,2)==0
%     disp('Order of polynomial must be odd. Stopping execution...');
%     return;
% end
Path = [];
WP   = [p0(1) p0(2); pt(1) pt(2)]
N    = length(WP)-1;

Path.NumSubpaths = 1;
Path.Order = ord;
Path.WP.x  = [p0(1); pt(1)];
Path.WP.y  = [p0(2); pt(2)];

%% Calculating subpaths 
A       = zeros(ord+1,ord+1);
coefs   = ones(1,ord+1);
A(1,1)  = coefs(1);
A(2,:)  = coefs;
c_der   = coefs;
for k=2:(ord+1)/2
    c_der       = fliplr(polyder(fliplr(c_der)));
    coefs       = [zeros(1,k-1) c_der];
    A(2*k-1,k)  = c_der(1);
    A(2*k,:)    = coefs;
end
Path.LinSys.A = A;

ax = zeros(ord+1,1);
bx = zeros(ord+1,1);
for j=1:N
    ax     = zeros(ord+1,1);
    bx     = zeros(ord+1,1);
    
    ax(1)  = WP(j,1);
    ax(2)  = WP(j+1,1);
    bx(1)  = WP(j,2);
    bx(2)  = WP(j+1,2);
    if ord>2
        ax(3) = lambda*t0(1)/norm(t0);
        bx(3) = lambda*t0(2)/norm(t0);
        ax(4) = lambda*(WP(j+1,1)-WP(j,1))/norm(pt-p0);
        bx(4) = lambda*(WP(j+1,2)-WP(j,2))/norm(pt-p0);
    end
    a_vec  = A\ax;
    b_vec  = A\bx;
    Path.LinSys.bx{j} = ax;
    Path.LinSys.by{j} = bx;
    Path.coeff.a{j}   = a_vec;
    Path.coeff.b{j}   = b_vec;
    a_vec = a_vec';
    b_vec = b_vec';
    for k=1:ord
        a_vec = fliplr(polyder(fliplr(a_vec)));
        Path.coeff.a_der{j}{k} = a_vec;
        b_vec = fliplr(polyder(fliplr(b_vec)));
        Path.coeff.b_der{j}{k} = b_vec;
    end
end


%% Plotting

if PlotHandle
    e = 0:.01:1;
    figure(PlotHandle); hold on;
    for j=1:N
        X = zeros(1,length(e));
        Y = zeros(1,length(e));
        for k=1:ord+1
            X = X + Path.coeff.a{j}(k)*e.^(k-1);
            Y = Y + Path.coeff.b{j}(k)*e.^(k-1);
        end

        for j=1:N+1
            plot(WP(j,2),WP(j,1),'r*','LineWidth',1.5); 
        end
        plot(Y,X,'b','LineWidth',1.25);
    end 
    hold off; xlabel('y-position'); ylabel('x-position');
%     MinX = min(WP(:,2)); MaxX = max(WP(:,2));
%     MinY = min(WP(:,1)); MaxY = max(WP(:,1));
    MinX = min(Y); MaxX = max(Y);
    MinY = min(X); MaxY = max(X);
    axis equal; axis([1.1*MinX-0.1*MaxX 1.1*MaxX-0.1*MinX 1.1*MinY-0.1*MaxY 1.1*MaxY-0.1*MinY]); 
end
