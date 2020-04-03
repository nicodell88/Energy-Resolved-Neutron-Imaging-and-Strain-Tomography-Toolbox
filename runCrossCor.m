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
opts.method     = 'crosscorr';
opts.order      = 3;            %polynomial order of Savitzky-Golay
opts.frame      = 9;            %frame-length for Savitzky-Golay (must be Odd)
opts.peakWindow = 20;           %window for peak fitting
opts.range      = [4.55 4.9];   %Wavelengths of interest

opts.p00 = [... %Pseudo-Voigt function params
    15.2407     %y0
    0.0         %xc
    1.8534      %log(A)
    -4.1025     %log(Wg)
    -3.5159     %log(Wl)
    0.0876];    %log(Mu)
%% Fit Edges
dd      = nan(15,1);
sigma   = nan(15,1);

%     [dd(i),sigma(i),stuff,things] = crossCorrMethod(Tr(:,i),d0_al,lambda,opts);
      [dd_cell,sigma_cell,stuff,things] = fitEdges({Tr.'},lambda,opts,d0_al.');

dd = dd_cell{1}
%% Plot Results
figure(3); clf;
errorbar(dist,dd,sigma)
xlabel('Z (distance from sample center) - [mm]')
ylabel('Shift in wavelength - [\AA]')
title('Al\{111\} Bragg Edge')
grid minor
