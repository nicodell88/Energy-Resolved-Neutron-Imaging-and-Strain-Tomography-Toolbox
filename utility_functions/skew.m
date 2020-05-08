function S = skew(u)
%skew    Obtain matrix operator for left action of cross product
%   S = skew(u)
%   takes column vector u of length 3 and returns a matrix S,
%   such that
%       S*v = cross(u,v)
%
%   See also cross.

S = [ ...
    0,-u(3),u(2); ...
    u(3),0,-u(1); ...
    -u(2),u(1),0 ...
    ];