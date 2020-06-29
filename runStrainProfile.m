%
% runStrainProfileExample.m
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


%% Options
profileOpts.ObFile     = '/Users/megatron/Dropbox/Strain Tomography 2020/LaTeX/ToolboxPaper/matlab/preprocessed/JPARC2018_OB_proj_003.mat';   %Open Beam File location
profileOpts.ProjFile   = '/Users/megatron/Dropbox/Strain Tomography 2020/LaTeX/ToolboxPaper/matlab/preprocessed/JPARC2018_proj_000.mat';    %Projection File location
% opts.ObFile     = '/Users/johannes/Dropbox/ToolboxPaper/matlab/preprocessed/JPARC2018_OB_proj_003.mat';   %Open Beam File location
% opts.ProjFile   = '/Users/johannes/Dropbox/ToolboxPaper/matlab/preprocessed/JPARC2018_proj_000.mat';    %Projection File location

OB      = load(profileOpts.ObFile);    %Load open beam
Proj    = load(profileOpts.ProjFile);  %Load projection
%% Options
profileOpts = ProfileOptions(Proj.tof,'sampleRange',[30 250]);
BraggOpts = BraggOptions(Proj.tof,'gp','Par',true);
%% Calculate Strain Profile
[StrainProfile,SigmaProfile,profileOpts]= makeStrainProfile(OB,Proj,profileOpts,BraggOpts);
%%
figure(2)
clf
errorbar(1e6*StrainProfile,1e6*SigmaProfile)
title('Strain Profile')
xlabel('X - [pixels]')
ylabel('Through Thickness Strain - [$\mu$-strain]') 

