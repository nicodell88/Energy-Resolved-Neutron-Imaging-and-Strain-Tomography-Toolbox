%
% runCompare.m 
%
% Copyright (C) 2020 The University of Newcastle, Australia
% Authors:
%   Nicholas O'Dell <Nicholas.Odell@newcastle.edu.au>
% Last modified: 10/01/2020
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
opts.range      = [0.018 0.019];    %Range for fitting 5 param method
opts.method     = 'gp';             %Fitting algorithm
opts.plot       = false;             %plot results along the way
%% Used for Attenuation method
opts.a00    = 0.5;          %Initial guess for a0
opts.b00    = 0.5;          %Initial guess for b0
opts.a_hkl0 = 0.5;          %Initial guess for a_hkl
opts.b_hkl0 = 0.5;          %Initial guess for b_hkl
%% Used for parametric methods
opts.t_hkl0     = 0.0187;      %Initial guess for edge location
opts.sigma0  	= 1e-5;        %Initial guess for gaussian broadening term
opts.tau0    	= 1e-5;        %Initial guess for exponential decay term
%% GP
opts.sig_f  = 1;            %Squared-Exponential Kernel Hyperparameter, output variance
opts.l      = 1e-4;         %Squared-Exponential Kernel Hyperparameter, lengthscale
opts.ns     = 3000;         %Number of MC samples used to estimate bragg-edge location and variance.
opts.n      = 2500;         %Number of points to sample the Bragg-Edge function.

testProj = {Tr{1}(55,:)};
opts.method     = 'attenuation';    %Fitting algorithm
[attenBragg,stdAtten,TrFit_cell1] = fitEdges(testProj,tof,opts);
opts.method     = '5param';         %Fitting algorithm
[fiveparamBragg,stdFive,TrFit_cell2] = fitEdges(testProj,tof,opts);
opts.method     = 'gp';             %Fitting algorithm
[gpBragg,stdGP,TrFit_cell3] = fitEdges(testProj,tof,opts);
%% Plot Results
figure(2); clf;
ax(1) = subplot(2,1,1);

plot(tof,TrFit_cell1{1});
hold on
plot(tof,TrFit_cell2{1});
plot(tof,TrFit_cell3{1});
plot(tof,testProj{1},'.');
xlabel('time-of-flight - [seconds]')
ylabel('Normalised Transmission')
legend('Attenuation Method','5 Parameter Method','GP Method','Data')
title('Edge Fitting Comparison')
ax(2) = subplot(2,1,2);
cla
x = linspace(tof(1),tof(end),2000);
plot(x,normpdf(x,attenBragg{1},stdAtten{1}))
hold on
plot(x,normpdf(x,fiveparamBragg{1},stdFive{1}))
plot(x,normpdf(x,gpBragg{1},stdGP{1}))
xlabel('time-of-flight - [seconds]')
ylabel('Density')
legend('Attenuation Method','5 Parameter Method','GP Method')
linkaxes(ax,'x');