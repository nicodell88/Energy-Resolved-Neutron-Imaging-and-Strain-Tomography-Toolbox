function [point, pos, isInOn] = intersectLineTri(line, triangle, varargin)
%INTERSECTLINETRIANGLE3D Intersection point of a 3D line and a 3D triangle
% lines is an 6xN with each coumn containing l0 and nhat

%% Default values
[r,N] = size(line);
if r ~= 6
    error('line should have 6 rows')
end

point = NaN(3,N);
pos = NaN(1,N);
isInOn = false(1,N);

tol = 1e-13;
if ~isempty(varargin)
    tol = varargin{1};
end


%% Process inputs

% triangle is given as a 3-by-3 array
t0  = triangle(:,1);
u   = triangle(:,2) - t0;
v   = triangle(:,3) - t0;

% triangle normal
n   = cross(u, v);
% test for degenerate case of flat triangle 
if norm(n) < tol
    return;
end

%% Compute intersection with plane
% % line direction vector
% dir = line(4:6,:);

% vector between triangle origin and line origin
w0 = line(1:3,:) - t0;

% compute projection of each vector on the plane normal
a = -dot(repmat(n,1,N), w0);
b = dot(repmat(n,1,N), line(4:6,:));

% test case of line parallel to the triangle and ignore these lines
inds = abs(b) > tol;        % maybe don't need this check as divide by 0 will give NaN??
NN = sum(inds);
if ~NN
   return; 
end
% compute intersection point of line with supporting plane
% If r < 0: point before ray
% If r > 1: point after edge
pos(inds) = a(inds)./ b(inds);

% coordinates of intersection point
point(:,inds) = line(1:3,inds) + repmat(pos(inds),3,1) .* line(4:6,inds);


%% test if intersection point is inside triangle

% normalize direction vectors of triangle edges
uu  = dot(u, u);
uv  = dot(u, v);
vv  = dot(v, v);

% coordinates of vector v in triangle basis
w   = point(:,inds) - t0;
wu  = dot(w, repmat(u,1,NN));
wv  = dot(w, repmat(v,1,NN));

% normalization constant
D = uv^2 - uu * vv;

% test first coordinate
s = (uv .* wv - vv .* wu) ./ D;
in1 = (s >= -eps) .* (s <=1.0+eps);

% test second coordinate, and third triangle edge
t = (uv .* wu - uu .* wv) ./ D;
in2 = (t >=-eps).*(s+t <= 1.0+eps);


isInOn(inds) = logical(in1.*in2);
point(:,~isInOn) = NaN;
pos(~isInOn) = NaN;

end

