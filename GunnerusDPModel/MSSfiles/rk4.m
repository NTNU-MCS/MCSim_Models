function xnext = rk4(f,x,u,h)
% RK4	Integrate a system of ordinary differential equations using
%	Runge-Kutta's 4th-order method.
%
% xnext = rk4(f,x,u,h)
%
% x     - x(k)
% u     - u(k)
% xnext - x(k+1)
% f     - function returning: dx/dt(k) = f(x(k),u(k))
% h     - step size
%
% Ex:   function dx = f(x,u), 
%          dx = sin(x)+u;
% ===>  xnext = rk4('f',x,u,h) 
%
% Author:   Thor I. Fossen
% Date:     14th June 2001
% Revisions: 

xo = x;
k1 = h*feval(f,xo,u);
x  = xo+0.5*k1;
k2 = h*feval(f,x,u);
x  = xo+0.5*k2;
k3 = h*feval(f,x,u);
x  = xo+k3;
k4 = h*feval(f,x,u);
xnext = xo + (k1+2*(k2+k3)+k4)/6;