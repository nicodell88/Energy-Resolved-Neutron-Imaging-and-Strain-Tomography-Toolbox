function [entry,exit,nhat,L,yInds,nsegs,Rejects] = find_intersects_2D(rBOo,rPOo,nhat)
%find_intersects_2D  Finds intersects between rays and the sample.
%[entry,exit,nhat,L,yInds,nsegs,Rejects] =
%find_intersects_2D(rBOo,rPOo,nhat)
%   Inputs:
%       - rBOo is a 2-by-N matrix, where each column is an x,y vertex
%           of the sample boundary. Closed polygons, annular shapes and
%           multi-body samples can be specified. E.g., a closed polygon, in 
%           this case a square, can be specified as:
%               rBOo = [0,0,1,1,0;0,1,1,0,0]; 
%           Annular or multi body shape boundaries are expressed by placing
%           a column of NaNs between the rows defining each shape. E.g., a
%           square within a square;
%               rBOo = [0,0,1,1,0,NaN,0.1,0.1,0.9,0.9,0.1;
%                       0,1,1,0,0,NaN,0.1,0.9,0.9,0.1,0.1];
%       - rPOo is a 2-by-M vector, each column defining the location
%           of each pixel associated with each measurement, in sample
%           coordinates, where M is the number of measurements.
%       - nhat is a 2-by-M matrix, where each column is a unit direction
%           vector defining the direction of the beam in sample
%           coordinates.
%
%   Outputs:
%       - entry is a G-by-H matrix defining where each beam enters the
%           sample, where H is the number of measurements and G is
%           two-times the maximum number of entrances into the sample, H is
%           the number of measurements. For convex polygons G=2, but for
%           non-convex polygons, annular shapes or multi-body samples G may
%           be greater. Where G is greater than 2 - columns associated with
%           measurements that correspond to rays that only enter the sample
%           once are padded with NaNs.
%       - exit is also a G-by-H matrix, defning where each beam leaves the
%           sample. See above.
%       - nhat, is simply returned.
%       - L is an H-by-1 vector, where each element is the total irradiated
%           length of each measurement.
%       - yInds is a vector of indicies indicating which measurements were
%           valid, i.e., if a ray does not intersect the sample at all, or
%           intersects the sample an odd number of times the measurement is
%           discarded.
%       - nSegs is a H-by-1 vector indicating how many disjoint line
%           segments make up the ray path.
%       - Rejects is a vector of indicies indicating which measurements
%           were deleted.
%
%   See also polyxpoly, GenerateLRT3D, GenerateLRT2D.

% Copyright (C) 2020 The University of Newcastle, Australia
% Authors:
%   Nicholas O'Dell <Nicholas.Odell@newcastle.edu.au>
%   Johannes Hendriks <Johannes.Hendriks@newcastle.edu.au>
% Last modified: 06/07/2020
% This program is licensed under GNU GPLv3, see LICENSE for more details.

n = length(nhat);

H = figure(1);
clf
plot(rPOo(1,:),rPOo(2,:))
hold on

X2 = rBOo(1,:);
Y2 = rBOo(2,:);

plot(X2,Y2)
quiver(rPOo(1,:),rPOo(2,:),nhat(1,:),nhat(2,:))
try
    intersections = cell(n,1);
    
    for i=1:n
        tmp = rPOo(:,i)+[-nhat(:,i),+nhat(:,i)];
        X1 = tmp(1,:);
        Y1 = tmp(2,:);
        
        [XI,YI] = polyxpoly(X1,Y1,X2,Y2);

        plot(XI,YI,'ro')
        drawnow
        
        intersections{i} = [XI(:).';YI(:).'];
    end
catch me
    close(H);
    rethrow(me);
end
close(H);
[~,numInts] = cellfun(@size,intersections);
maxInts = max(numInts);

entry = nan(ceil(maxInts)/2*2,n);
exit = nan(ceil(maxInts)/2,n);
L = zeros(n,1);
for k =1:n
    N = size(intersections{k},2)/2*2;
    entry(1:2:N,k) = intersections{k}(1,1:2:end);%Xentry
    entry(2:2:N,k) = intersections{k}(2,1:2:end);%Yentry
    
    exit(1:2:N,k) = intersections{k}(1,2:2:end);%Xexit
    exit(2:2:N,k) = intersections{k}(2,2:2:end);%Xexit
    
    L(k) = 0;
    for i=1:N/2
        L(k) = L(k) + norm(entry(1:2+(i-1),k)-exit(1:2+(i-1),k));
    end

end
del = logical(mod(numInts,2));
yInds = find(~del);
Rejects = find(del);
nsegs = round(numInts/2);

nsegs(del)=[];
entry(:,del) = [];
exit(:,del) = [];


end