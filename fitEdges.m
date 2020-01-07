function [d_cell,std_cell] = fitEdges(Tr,tof,opts)
%FITEDGES Fits Bragg edges to input data.
%   Inputs:
%       - Tr is a cell-array of normalised transmission intensity data. Each cell
%       is a projection.
%       - tof is an array of wave-lengths of time-of-flight, whos length matches
%       the width of each array in Tr.
%       - options is a structure containing
%           opts.method      :  a char array defining the fitting method 
%                               ('Santisteban2001', 'Tremsin2011', 'Hendriks2020')
%           opts.startRange  :  a 2 element vector containing the start 
%                               and end range for fitting the left side of the 
%                               Bragg-Edge [start end]
%           opts.endRange    :  a 2 element vector containing the start 
%                               and end range for fitting the right side of 
%                               the Bragg-Edge [start end]
%   Outputs:
%       - d_cell is a cell array containing the bragg edge location for
%       each projection.
%       - std_cell is a cell array containing corresponding standard
%       deviation estimates for each result in d_cell.
%
% Copyright (C) 2020  University of Newcastle, Australia
% Authors:
%   Nicholas O'Dell <Nicholas.Odell@newcastle.edu.au>
% Last modified: 07/01/2020
% This program is licensed under GNU GPLv3, see LICENSE for more details.

% least squares fitting options
optionsFit = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');
optionsFit.Algorithm = 'Levenberg-Marquardt';
optionsFit.Jacobian = 'off';
optionsFit.Display = 'off';

%TODO write function
end

