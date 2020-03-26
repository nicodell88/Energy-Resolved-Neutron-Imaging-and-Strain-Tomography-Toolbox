function [x] = unscaleInput(xscaled,xrange)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

x = xscaled * (xrange(2)-xrange(1)) + xrange(1);

end

