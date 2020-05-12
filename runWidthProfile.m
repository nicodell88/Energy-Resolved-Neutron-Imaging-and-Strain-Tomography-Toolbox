% runWidthProfileExample.m 
%
% Copyright (C) 2020 The University of Newcastle, Australia
% Authors:
%   Nicholas O'Dell <Nicholas.Odell@newcastle.edu.au>
%   Johannes Hendriks <Johannes.Hendriks@newcastle.edu.au>
% Last modified: 12/05/2020
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


%% Options
% opts.ObFile     = '/Users/megatron/Dropbox/Strain Tomography 2020/J-PARC 2020 Working Folder/Data/OpenBeam/OpenBeamImg__proj001_idx.mat';   %Open Beam File location
% opts.ProjFile   = '/Users/megatron/Dropbox/Strain Tomography 2020/J-PARC 2020 Working Folder/Data/SteelCube/SteelCube__proj001_idx.mat';    %Projection File location
opts.ObFile     = '/Users/johannes/Dropbox/ToolboxPaper/matlab/preprocessed/JPARC2018_OB_proj_003.mat';   %Open Beam File location
opts.ProjFile   = '/Users/johannes/Dropbox/ToolboxPaper/matlab/preprocessed/JPARC2018_proj_000.mat';    %Projection File location


OB      = load(opts.ObFile);    %Load open beam
Proj    = load(opts.ProjFile);  %Load projection
%% Options
% opts.supplyMask = 'supply';
% opts.supplyMask = 'gui';
opts.supplyMask = 'auto';

opts.rangeRight = [0.0186 0.0188];  %Range for averaging on right of edge for guessing mask
opts.rangeLeft  = [0.0175 0.01835];  %Range for averaging on left of edge for guessing mask
% opts.d0         =  0.018821361897845;           %Unstrained tof/wl
opts.d0         =  0.01845;           %Unstrained tof/wl
%% PROFILE OPTIONS
opts.pixelCentre = 120;
opts.nWidth     = 3;    %
opts.nPix       = 80;   %averages over (2*nPix+1)-by-(2*nWidth+1)
opts.nRes       = 2;
%%
opts.Thresh     = 0.1;             %More than 90% of a macro-pixel must be within the mask before a bragg edge will be fit to it
opts.maskThresh = 0.02;             %Threshold edge height for mask


opts.AveDir      = 'horz'; %'horz'

% Insert Bragg Edge Fitting Options here
opts.BraggOpts.startRange = [0.0182 0.01835];  %Fitting left side of edge
opts.BraggOpts.endRange   = [0.0186 0.0188];   %Fitting right side of edge

opts.BraggOpts.GPscheme     = 'hilbertspace';             %Fitting algorithm
% opts.BraggOpts.method = 'gp';
opts.BraggOpts.method = 'attenuation';
opts.BraggOpts.covfunc = 'M52';        %Sets the GP covariance function (currently only SE implemented)
opts.BraggOpts.optimiseHP = 'all';

opts.BraggOpts.sigma0     = 0.01;        %Initial guess for gaussian broadening term
opts.BraggOpts.tau0       = 0.01;        %Initial guess for exponential decay term
opts.t_hkl0     = 0.01845;     %Initial guess for edge location

opts.BraggOpts.plot   = false;             %plot results along the way
opts.BraggOpts.a00    = 0.5;          %Initial guess for a0
opts.BraggOpts.b00    = 0.5;          %Initial guess for b0
opts.BraggOpts.a_hkl0 = 0.5;          %Initial guess for a_hkl
opts.BraggOpts.b_hkl0 = 0.5;          %Initial guess for b_hkl
opts.BraggOpts.sig_f  = 1;            %Squared-Exponential Kernel Hyperparameter, output variance
opts.BraggOpts.l      = 0.0445;         %Squared-Exponential Kernel Hyperparameter, lengthscale
opts.BraggOpts.ns     = 500;         %Number of MC samples used to estimate bragg-edge location and variance.
opts.BraggOpts.n      = 2500;         %Number of points to sample the Bragg-Edge function.




% Calculate Strain Profile
opts.BraggOpts.Par = true;
opts.figNum = 1;
% tic
[strainProfile,sigmaProfile,opts,widthProfile,widthStdProfile]= makeStrainProfile(OB,Proj,opts);
%%
inds = 40:(length(strainProfile)/2-40);
figure(1)
clf
subplot 211
errorbar(inds,strainProfile(inds),sigmaProfile(inds))
title('Strain profile')
subplot 212
errorbar(inds,widthProfile(inds),widthStdProfile(inds))
title('Width profile')
% ylim([0 0.1])
% xlim([min(inds) max(inds)])
% toc
%% Plot Strain Image
% plotStrainImage(StrainProfile,SigmaProfile,opts)
%% Save Figure
% msg = sprintf('strainImage_%s',datestr(now,'yy_mm_dd_hh_MM_ss'));
% saveas(gcf,msg,'png')
