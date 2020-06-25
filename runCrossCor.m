%
% runCrossCor.m
%
% Copyright (C) 2020 The University of Newcastle, Australia
% Authors:
%   Nicholas O'Dell <Nicholas.Odell@newcastle.edu.au>
% Last modified: 14/01/2020
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
load    corr_paper_data.mat

%% Fit Bragg-Edge
opts = BraggOptions(lambda,'crosscorr','frame',35,'peakwindow',16,'Wg',1000,'Wl',10,'Mu',0.1);

%% Fit Edges
[dd_cell,sigma_cell,~,fitinfo] = fitEdges({Tr.'},lambda,opts,d0_al.');

%% Plot Results
figure(3); clf;
errorbar(dist,dd_cell{1}.',sigma_cell{1})
xlabel('Z (distance from sample center) - [mm]')
ylabel('Shift in wavelength - [\AA]')
title('Al\{111\} Bragg Edge')
grid minor
