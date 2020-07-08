function [processed] = ProcessMat(input,varargin)
% Processes a .mat file (I0 or I) and provides back an intensity profile.
% One function for both to ensure consistency.
%
% Code Structure: 
% - Convert TOF to wavelength.
% - Sum over all columns to form a 2D detector.
% - ./ by ntrigs to account for non-uniform sample times.
% - trims pixels outside the sample range

% Copyright (C) 2020 The University of Newcastle, Australia
% Authors:
%   Nicholas O'Dell <Nicholas.Odell@newcastle.edu.au>
%   Johannes Hendriks <Johannes.Hendriks@newcastle.edu.au>
% Last modified: 06/07/2020
% This program is licensed under GNU GPLv3, see LICENSE for more details.

%%
if nargin == 1
% From Anton
trig_delay = 1.243e-5; % [S]
% From J-PARC calibration sample
source_dist = 17.7971; % [m] 
elseif nargin == 3
        trig_delay = varargin{1};
        validateattributes(trig_delay,{'numeric'},{'scalar','positive'})
        source_dist = varargin{2};
        validateattributes(source_dist,{'numeric'},{'scalar','positive'})
else
    error('Must provide either both trigger delay and source distance. Or neither.')
end

%% - Convert TOF to wavelength.
% Fundamental Parameters of the Universe
plank = 6.62607004e-34; % [m^2 kg /s]
neutron_mass = 1.6749274e-27; % [kg]

% Calculate the wavelengths
lambda = plank/neutron_mass/source_dist*(input.tof-trig_delay)*1e10;

%% dot-divide intensities by ntrigs to account for non-uniform sample times
% For example, open-beam (or d0) might be much longer than a single projection.
% Dividing by ntrigs gives back a relative intensity that can be compared.
ntrigs_rep = reshape(input.ntrigs,1,1,numel(input.ntrigs)).*ones(size(input.im_stack,1),size(input.im_stack,2),numel(input.ntrigs));

%% Calculate outputs
processed.stack = input.im_stack./ntrigs_rep;
processed.lambda = lambda;

% disp('done.')
end