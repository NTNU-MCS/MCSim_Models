function S = Smtrx(r)
% S = SMTRX(r) computes the 3x3 vector cross product matrix S=-S'
% such that rxt = S(r)t is true for all 3x1 vectors r and t
%
% Author:   Thor I. Fossen
% Date:     14th June 2001
% Revisions: 
 
S = [    0  -r(3)   r(2)
      r(3)     0   -r(1)
     -r(2)   r(1)     0 ];
 