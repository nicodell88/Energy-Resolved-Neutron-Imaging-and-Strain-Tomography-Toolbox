function [xscaled,xrange] = scaleInput(x)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
xmax = max(x);
xmin = min(x);
xrange = [xmin,xmax];

if abs(xmax - xmin) < sqrt(eps)
    warning("Not a valid range")
    xscaled = x;
end

xscaled = (x - xmin)/(xmax-xmin);

    
% xscaled = x;
% xrange = nan;

end

