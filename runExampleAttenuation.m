%
% runExampleAttenuation.m 
%
% Copyright (C) 2020 The University of Newcastle, Australia
% Authors:
%   Nicholas O'Dell <Nicholas.Odell@newcastle.edu.au>
% Last modified: 20/04/2020
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.

restoredefaultpath
clc
clear
close all
%% House Keeping
tickFontSize = 12;
set(0,'defaultAxesFontsize',tickFontSize)
set(0,'defaultAxesLabelFontsize',18/tickFontSize)
set(0,'defaultAxesTitleFontsize',25/tickFontSize)
set(0,'defaultTextInterpreter','latex')
set(0,'defaultLegendInterpreter', 'latex')
set(0,'defaultAxesTickLabelInterpreter', 'latex')
set(0, 'defaultLineLinewidth',2)
set(0,'defaultFigureColor','w')
%% Load Data
addpath ./data
addpath ./utility_functions
load CubeProj.mat
%% Plot Data
figure(1); clf;
plot(tof,Tr{3}(100,:))%arbitrary
xlabel('Time-Of-Flight - [seconds]')
ylabel('Normalised Transmission Intensity - [arbitrary units]')
grid minor
%% Fit Bragg-Edge
opts.startRange = [tof(1) tof(150)];    %Fitting left side of edge
opts.endRange   = [tof(371) tof(end)];  %Fitting right side of edge
opts.method     = 'attenuation';    %Fitting algorithm
opts.plot       = true;                 %plot results along the way

opts.a00        = 0.5;          %Initial guess for a0
opts.b00        = 0.5;          %Initial guess for b0
opts.a_hkl0     = 0.5;          %Initial guess for a_hkl
opts.b_hkl0     = 0.5;          %Initial guess for b_hkl
opts.t_hkl0     = tof(250);     %Initial guess for edge location
opts.sigma0     = log(4.4e-5);        %Initial guess for gaussian broadening term
opts.tau0       = 0.008;        %Initial guess for exponential decay term

% [d_cell,std_cell,TrFit_cell,fitInfo_cell] = fitEdges(Tr,tof,opts);

opts.Par = false;
% [d_cell,std_cell,TrFit_cell,fitInfo_cell] = fitEdges(Tr(1),tof,opts);
[d_cell,std_cell,TrFit_cell,fitInfo_cell] = fitEdges(Tr,tof,opts);
