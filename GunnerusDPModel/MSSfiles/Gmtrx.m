function G = Gmtrx(nabla,A_wp,GMT,GML,r_g)
% G = GMTRX(nabla,A_wp,GMT,GML,r_g) computes the 6x6 system spring stiffness matrix G
% about an arbitrarily point O for a floating vessel (small roll and pitch angles).
% For submerged vessels, see gvect.m
% 
% Inputs:  nabla: deplasement
%          Awp: water plane area
%          GMT, GML: transverse/longitudinal metacentric heights
%          r_g = [x_g y_g z_g]': location of CG with respect to O
%
% Author:     Thor I. Fossen
% Date:       14th June 2001
% Revisions:  26th June 2002,  variable Awp was replaced with A_wp 
%                              one zero in G_CG was removed

rho = 1025;  % density of water
g   = 9.81;	 % acceleration of gravity

Zz     = -rho*g*A_wp;
Kphi   = -rho*g*nabla*GMT;
Mtheta = -rho*g*nabla*GML;
G_CG   = diag([0 0 -Zz -Kphi -Mtheta 0]);
G = Hmtrx(r_g)' * G_CG * Hmtrx(r_g);