function [entry,exit,nhat,L,yInds,nsegs,rejects] = find_intersects_3D(shape,lines)
% function [entry,exit,nhat,L,yInds,nsegs,rejects] = find_intersects_3D_exp(shape,lines)
%[entry,exit,nhat,L,yInds,nsegs,rejects] = find_intersects_3D(shape,lines)
%
%Inputs:
%   - shape is a structure with
%       N: normals to each face, N(i,:) is the normal to the ith face.
%       C: centers of each face, C(i,:) is the center of the ith face.
%       F: Indices of vertices corresponding to each face, F(i,:) are the
%       indices of the vertices corresponding to face i.
%       V: vertices of the shape, V(k,:) = [x,y,z] of the kth vertices.
%   - lines is a 2D array containing the starting and direction of each line,
%       line(:,j) = [x,y,z,nhat_x,nhat_y,nhat_z]' is the start point and
%       direction of the jth line.
%Outputs:
%   - entry is a G-by-H matrix defining where each beam enters the
%       sample, where H is the number of measurements and G is
%       two-times the maximum number of entrances into the sample, H is
%       the number of measurements. For convex polygons G=2, but for
%       non-convex polygons, annular shapes or multi-body samples G may
%       be greater. Where G is greater than 3 - columns associated with
%       measurements that correspond to rays that only enter the sample
%       once are padded with NaNs.
%   - exit is also a G-by-H matrix, defning where each beam leaves the
%       sample. See above.
%   - nhat, is simply returned.
%   - L is an H-by-1 vector, where each element is the total irradiated
%       length of each measurement.
%   - yInds is a vector of indicies indicating which measurements were
%       valid, i.e., if a ray does not intersect the sample at all, or
%       intersects the sample an odd number of times the measurement is
%       discarded.
%   - nSegs is a H-by-1 vector indicating how many disjoint line
%       segments make up the ray path.
%   - Rejects is a vector of indicies indicating which measurements
%       were deleted.
%
%See also downSample3D_LRT.

% Copyright (C) 2020 The University of Newcastle, Australia
% Authors:
%   Johannes Hendriks <Johannes.Hendriks@newcastle.edu.au>
%   Nicholas O'Dell <Nicholas.Odell@newcastle.edu.au>
% Last modified: 16/07/2020
% This program is licensed under GNU GPLv3, see LICENSE for more details.

[nn,N] = size(lines);
assert(nn==6,'line(:,k) must constain [x,y,z,nhat_x,nhat_y,nhat_z] transposed');

[num_faces,~] = size(shape.F);

Point = NaN(3,num_faces*N);
Pos = NaN(1,N*num_faces);
IsInOn = NaN(1,N*num_faces);
rayInd = repmat(1:N,1,num_faces);
for i = 1:num_faces
    tri = shape.V(shape.F(i,:),:)';
    [Point(:,N*i-(N-1):N*i), Pos(N*i-(N-1):N*i), IsInOn(N*i-(N-1):N*i)] = intersectLineTri(lines, tri);
end
% remove NaNs
Point(:,~IsInOn) = [];
Pos(~IsInOn) =[];
rayInd(~IsInOn) =[];

% remove duplicates
[C,IA,IC] = unique(Point','rows','stable');
%
Point = C';
Pos = Pos(IA);
rayInd = rayInd(IA);

% sort by distance along ray
[s, sInd] = sort(Pos);
Point = Point(:,sInd);
rayInd = rayInd(sInd);

% sort by ray ind, this must be done after sorting by distance
[rayInd,rInd] = sort(rayInd);
Point = Point(:,rInd);
s = s(rInd);

% % check number of intersections per ray
intsPerRay = hist(rayInd,1:N);
% uneven number of intersects are caught later

maxIntsPerRay = max(intsPerRay);
if max(intsPerRay) <=2 % this stuff only works if all rays have 0 or 2 intersections
    magicInd = rayInd*max(intsPerRay)-(max(intsPerRay)-1) + [0 ~diff(rayInd)];
else
    magicInd = NaN(N*max(intsPerRay),1);
    prior_ints = 0;
    start_ind = 1;
    for ii = 1:N
        magicInd(start_ind:start_ind+intsPerRay(ii)-1) = prior_ints+1:prior_ints+intsPerRay(ii);
        prior_ints = prior_ints+maxIntsPerRay;
        start_ind = start_ind + intsPerRay(ii);
    end
    
end
magicInd = magicInd(~isnan(magicInd));

tmpx = NaN(max(intsPerRay),N);
tmpy = NaN(max(intsPerRay),N);
tmpz = NaN(max(intsPerRay),N);

tmpx(magicInd) = Point(1,:)';
tmpy(magicInd) = Point(2,:)';
tmpz(magicInd) = Point(3,:)';

entry = NaN(ceil(max(intsPerRay)/2)*3,N);
exit = NaN(ceil(max(intsPerRay)/2)*3,N);

entry(1:3:end,:) = tmpx(1:2:end,:);
entry(2:3:end,:) = tmpy(1:2:end,:);
entry(3:3:end,:) = tmpz(1:2:end,:);

exit(1:3:end,:) = tmpx(2:2:end,:);
exit(2:3:end,:) = tmpy(2:2:end,:);
exit(3:3:end,:) = tmpz(2:2:end,:);

del = logical( ~intsPerRay + mod(intsPerRay,2));
yInds = find(~del);
entry(:,del) = [];
exit(:,del) = [];
rejects = find(del);
nsegs = round(intsPerRay/2); nsegs(del) =[];


L = zeros(N,1);
for k=1:N
    
    for i = 1:nsegs(k)
        idx = (1:3)+(i-1)*3;
        L(k) = L(k) + norm(exit(idx,k) - entry(idx,k),2);
    end
end
nhat = lines(4:6,:);

end



