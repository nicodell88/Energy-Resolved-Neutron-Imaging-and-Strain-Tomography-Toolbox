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
opts.ProjFile   = '/Users/megatron/Dropbox/Strain Tomography 2020/J-PARC 2020 Working Folder/Data/SteelCube/SteelCube__proj001_idx.mat';    %Projection File location

OB      = load(opts.ObFile);    %Load open beam
Proj    = load(opts.ProjFile);  %Load projection
%%
% opts.supplyMask = 'supply';
opts.supplyMask = 'gui';
% opts.supplyMask = 'auto';

opts.rangeRight = [Proj.tof(end-50) Proj.tof(end)];  %Range for averaging on right of edge for guessing mask
opts.rangeLeft  = [Proj.tof(1) Proj.tof(50)];  %Range for averaging on left of edge for guessing mask
opts.d0         =  0.018821361897845;           %Unstrained tof/wl
opts.nPix       = 18;               %number of pixels to average over (npix-by-npix)
opts.nRes       = 18;               %number of pixels to step by before calculating next edge
opts.Thresh     = 0.01;             %More than 99% of a macro-pixel must be within the mask before a bragg edge will be fit to it
opts.maskThresh = 0.02;             %Threshold edge height for mask

opts.BraggOpts = BraggOptions(Proj.tof,'gp','Par',true);

[StrainImage,SigmaImage,opts]= makeStrainImage(OB,Proj,opts);

%% Plot Strain Image
plotStrainImage(StrainImage-nanmean(StrainImage,'all'),SigmaImage,opts)
%% Save Figure
msg = sprintf('strainImage_%s',datestr(now,'yy_mm_dd_hh_MM_ss'));
saveas(gcf,msg,'png')






