function [rPDd] = pix2vec(u,v)
%[rPDd] = pix2vec(u,v) Converts a list of pixels defined by u and v to the
%vector rPDd, returned as a 3-by-N vector where N is the number of pixels.
%
% Copyright (C) 2020 The University of Newcastle, Australia
% Authors:
%   Nicholas O'Dell <Nicholas.Odell@newcastle.edu.au>
% Last modified: 25/01/2020
% This program is licensed under GNU GPLv3, see LICENSE for more details.

assert(all(u<513 & u>0,'all'),'Pixel coordinates should be in the range [1,512]')
assert(all(v<513 & v>0,'all'),'Pixel coordinates should be in the range [1,512]')
assert(numel(u) == numel(v));

u = u(:).';
v = v(:).';

n = numel(u);

rPDd = [...
    zeros(1,n);
    55e-6*(v-513/2);
    55e-6*(u-513/2)];


end

