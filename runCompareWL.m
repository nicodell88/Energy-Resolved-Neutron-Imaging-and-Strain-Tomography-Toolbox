%
% runCompareWL.m 
%
% Copyright (C) 2020 The University of Newcastle, Australia
% Authors:
%   Nicholas O'Dell <Nicholas.Odell@newcastle.edu.au>
%   Johannes Hendriks <johannes.hendriks@newcastle.edu.au>
% Last modified: 18/03/2020
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
plot(wl,Tr{1}(20,:))%arbitrary
xlabel('Wave Length - [\AA]')
ylabel('Normalised Transmission Intensity - [arbitrary units]')
grid minor
%% Fit Bragg-Edge
attenOpts       = BraggOptions(wl,'attenuation');
fiveParamOpts   = BraggOptions(wl,'5param','range',[4.1 4.2]);
gpOpts          = BraggOptions(wl,[],'l',0.1);

testProj = {Tr{1}(200,:)};
[attenBragg,stdAtten,TrFit_cell1] = fitEdges(testProj,wl,attenOpts);
[fiveparamBragg,stdFive,TrFit_cell2] = fitEdges(testProj,wl,fiveParamOpts);
[gpBragg,stdGP,TrFit_cell3] = fitEdges(testProj,wl,gpOpts);

%% Plot Results
figure(2); clf;
ax(1) = subplot(2,1,1);

plot(wl,TrFit_cell1{1});
hold on
plot(wl,TrFit_cell2{1});
plot(wl,TrFit_cell3{1});
plot(wl,testProj{1},'.');
xlabel('Wave Length - [\AA]')
ylabel('Normalised Transmission')
legend('Attenuation Method','5 Parameter Method','GP Method','Data')
title('Edge Fitting Comparison')
ax(2) = subplot(2,1,2);
cla
x = linspace(wl(1),wl(end),2000);
plot(x,normpdf(x,attenBragg{1},stdAtten{1}))
hold on
plot(x,normpdf(x,fiveparamBragg{1},stdFive{1}))
plot(x,normpdf(x,gpBragg{1},stdGP{1}))
xlabel('Wave Length - [\AA]')
ylabel('Probability Density')
legend('Attenuation Method','5 Parameter Method','GP Method')
linkaxes(ax,'x');