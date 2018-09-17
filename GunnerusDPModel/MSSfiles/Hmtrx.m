function H = Hmtrx(r)
% H = HMTRX(r) computes the 6x6 system transformation matrix
%
% Author:   Thor I. Fossen
% Date:     14th June 2001
% Revisions: 
 
S = Smtrx(r);
H = [eye(3)     S'
     zeros(3,3) eye(3) ];
