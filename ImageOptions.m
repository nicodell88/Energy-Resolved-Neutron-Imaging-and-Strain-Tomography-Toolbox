function [opts] = ImageOptions(tof,varargin)
%ImageOptions(tof,...) generates an options structure which contains the
%neccecary options for strain imaging
%   opts = ImageOptions(tof) specifies the options for downsampling the
%   data using 18x18 macro pixels and will automatically generate the mask.   
%
%   Namve-value pair arguments can be used to change any of the parameters
%   from their default value. see makeStrainImage to see what they are.
%
%   See also makeStrainImage, BraggOptions.
%

% Copyright (C) 2020 The University of Newcastle, Australia
% Authors:
%   Nicholas O'Dell <Nicholas.Odell@newcastle.edu.au>
% Last modified: 24/06/2020
% This program is licensed under GNU GPLv3, see LICENSE for more details.

validateattributes((tof),{'numeric'},{'increasing'})
p = inputParser;


%% Add options
%% supplyMask
default = 'auto';
expectedMethods = {'auto','gui','supply'};
addOptional(p,'supplyMask',default, ...
    @(x) any(validatestring(x,expectedMethods)));
%% rangeLeft
default = [tof(1) tof(30)];
addParameter(p,'rangeLeft',default, ...
    @(x) validateattributes((x),{'numeric'},{'vector','numel',2,'increasing','>=',tof(1),'<=',tof(end)}));
%% endRange
default = [tof(end-30) tof(end)];
addParameter(p,'rangeRight',default, ...
    @(x) validateattributes((x),{'numeric'},{'vector','numel',2,'increasing','>=',tof(1),'<=',tof(end)}));
%% d0
default = tof(round(length(tof)/2));
addParameter(p,'d0',default, ...
    @(x) validateattributes((x),{'numeric'},{'scalar','>=',tof(1),'<=',tof(end)}));
%% nPix
default = 18;
addParameter(p,'nPix',default, ...
    @(x) validateattributes((x),{'numeric'},{'scalar','positive','integer'}));
%% nRes
default = 18;
addParameter(p,'nRes',default, ...
    @(x) validateattributes((x),{'numeric'},{'scalar','positive','integer'}));
%% Thresh
default = 0.01;
addParameter(p,'Thresh',default, ...
    @(x) validateattributes((x),{'numeric'},{'scalar','positive','<=',1}));
%% maskThresh
default = 0.02;
addParameter(p,'maskThresh',default, ...
    @(x) validateattributes((x),{'numeric'},{'scalar'}));
%%
parse(p,varargin{:});
opts = p.Results;

end

