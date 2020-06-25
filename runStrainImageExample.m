%
% runStrainImageExample.m 
%
% Copyright (C) 2020 The University of Newcastle, Australia
% Authors:
%   Nicholas O'Dell <Nicholas.Odell@newcastle.edu.au>
% Last modified: 25/06/2020
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
ObFile     = '/Users/megatron/Dropbox/Strain Tomography 2020/J-PARC 2020 Working Folder/Data/OpenBeam/OpenBeamImg__proj001_idx.mat';   %Open Beam File location
ProjFile   = '/Users/megatron/Dropbox/Strain Tomography 2020/J-PARC 2020 Working Folder/Data/SteelCube/SteelCube__proj001_idx.mat';    %Projection File location

OB      = load(ObFile);    %Load open beam
Proj    = load(ProjFile);  %Load projection

%% Options
BraggOpts = BraggOptions(Proj.tof,'gp','Par',true);
ImageOpts = ImageOptions(Proj.tof);

%% Generate strain image
[StrainImage,SigmaImage,opts]= makeStrainImage(OB,Proj,ImageOpts,BraggOpts);

%% Plot Strain Image
plotStrainImage(StrainImage-nanmean(StrainImage,'all'),SigmaImage,ImageOpts)
%% Save Figure
msg = sprintf('strainImage_%s',datestr(now,'yy_mm_dd_hh_MM_ss'));
saveas(gcf,msg,'png')






