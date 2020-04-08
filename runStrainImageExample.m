%
% runStrainImageExample.m 
%
% Copyright (C) 2020 The University of Newcastle, Australia
% Authors:
%   Nicholas O'Dell <Nicholas.Odell@newcastle.edu.au>
% Last modified: 05/03/2020
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
opts.ObFile     = '/Users/megatron/Dropbox/Strain Tomography 2020/J-PARC 2020 Working Folder/Data/OpenBeam/OpenBeamImg__proj001_idx.mat';   %Open Beam File location
opts.ProjFile   = '/Users/megatron/Dropbox/Strain Tomography 2020/J-PARC 2020 Working Folder/Data/SteelCube/SteelCube__proj050_idx.mat';    %Projection File location

% opts.ObFile = '/Users/megatron/Dropbox/Strain Tomography 2020/J-PARC 2019 Working Folder/data/open_beam/open_beam_11ksec.mat';
% opts.ProjFile = '/Users/megatron/Dropbox/Strain Tomography 2020/J-PARC 2019 Working Folder/data/cube_plug/Cube1__proj54_idx.mat';

% opts.ProjFile   = './Data/SteelCube/SteelCube__proj051_idx.mat';    %Projection File location
% opts.ProjFile = '/Users/megatron/Dropbox/Strain Tomography 2020/J-PARC 2020 Working Folder/Data/SteelCube/SteelCube__proj060_idx.mat'

OB      = load(opts.ObFile);    %Load open beam
Proj    = load(opts.ProjFile);  %Load projection
%%
% opts.supplyMask = 'supply';
% opts.supplyMask = 'gui';
opts.supplyMask = 'auto';

opts.rangeRight = [0.0190 0.0198];  %Range for averaging on right of edge for guessing mask
opts.rangeLeft  = [0.0175 0.0185];  %Range for averaging on left of edge for guessing mask
opts.d0         =  0.018821361897845;           %Unstrained tof/wl
% opts.d0         =  0.018821361897845 + 3e-6;           %Unstrained tof/wl
opts.nPix       = 18;               %Macro-pixel size
opts.Thresh     = 0.2;             %More than 80% of a macro-pixel must be within the mask before a bragg edge will be fit to it
opts.maskThresh = 0.02;             %Threshold edge height for mask
% Insert Bragg Edge Fitting Options here
opts.BraggOpts.startRange = [0.0175 0.0180];  %Fitting left side of edge
opts.BraggOpts.endRange   = [0.019 0.0195];   %Fitting right side of edge
% opts.BraggOpts.method     = 'attenuation';             %Fitting algorithm
opts.BraggOpts.GPscheme     = 'hilbertspace';             %Fitting algorithm
%  
opts.BraggOpts.method = 'gp';
% opts.BraggOpts.method = 'attenuation';
opts.BraggOpts.covfunc = 'M52';        %Sets the GP covariance function (currently only SE implemented)
opts.BraggOpts.optimiseHP = 'none';



opts.BraggOpts.plot   = false;             %plot results along the way
opts.BraggOpts.a00    = 0.5;          %Initial guess for a0
opts.BraggOpts.b00    = 0.5;          %Initial guess for b0
opts.BraggOpts.a_hkl0 = 0.5;          %Initial guess for a_hkl
opts.BraggOpts.b_hkl0 = 0.5;          %Initial guess for b_hkl
opts.BraggOpts.sig_f  = 1;            %Squared-Exponential Kernel Hyperparameter, output variance
opts.BraggOpts.l      = 0.04;         %Squared-Exponential Kernel Hyperparameter, lengthscale
opts.BraggOpts.ns     = 3000;         %Number of MC samples used to estimate bragg-edge location and variance.
opts.BraggOpts.n      = 2500;         %Number of points to sample the Bragg-Edge function.

%% Calculate Strain Image
opts.figNum = 1;
[StrainImage,SigmaImage,opts]= makeStrainImage(OB,Proj,opts);

 %Create temporary variables to mess with
StrainImagePlot = StrainImage;
% pcolor(StrainImage)
% shading flat

%% Attempting to manually removing outliers
% thresh = 1e-2
% idx = StrainImage > thresh | StrainImage < -thresh;
% StrainImagePlot = StrainImage;
% StrainImagePlot(idx) = nan;
% StrainImagePlot = StrainImagePlot - nanmean(StrainImagePlot,'all')
%% Plot Strain Image
plotStrainImage(StrainImagePlot,SigmaImage,opts)
%% Save Figure
msg = sprintf('strainImage_%s',datestr(now,'yy_mm_dd_hh_MM_ss'));
% saveas(gcf,msg,'png')
%% Here is where I was returning d0 in the StrainImage output an d in the sigma image output
figure(2)
clf
subplot(1,2,1)
pcolor(opts.d0Image)
shading flat
title('$\lambda_0$')
colorbar
subplot(1,2,2)
pcolor(opts.dImage)
shading flat
title('$\lambda$')
colorbar
%% Here is where I was returning d0 in the StrainImage output an d in the sigma image output
% figure(3)
% clf
% A = (SigmaImage-StrainImage)./StrainImage
% pcolor(A)
% shading flat
% colorbar
% caxis([0 8e-3])





