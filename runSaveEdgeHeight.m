%
% runSaveEdgeHeight.m
%
% Copyright (C) 2020 The University of Newcastle, Australia
% Authors:
%   Nicholas O'Dell <Nicholas.Odell@newcastle.edu.au>
% Last modified: 21/05/2020
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
%%
addpath ../utility_functions
%%
opts.rangeLeft  = [0.01825 0.01846];        %range for averaging on the left of the edge (same units as projX/tof)
opts.rangeRight = [0.01846 0.01912];        %range for averaging on the right of the edge (same units as projX/tof)
opts.fmt        = 'Cube1__proj%d_idx.mat';  %data file nmae spec
opts.proj_idx   = [0:69];                   %projections to be used for sample alignment
projPath = uigetdir('/Users/megatron/Dropbox/Strain Tomography 2020/J-PARC 2019 Working Folder/data/cube_plug');
opts.OB  = '/Users/megatron/Dropbox/Strain Tomography 2020/J-PARC 2019 Working Folder/data/open_beam/open_beam_11ksec.mat';
tic
edges = getEdgeHeight(projPath,opts);
toc
%%
figure(1)
clf
xlabel('hi')
for i = 1:size(edges,3)
    imagesc(edges(:,:,i))
    xlabel('X - [pixels]')
    ylabel('Y - [pixels]')
    title('Edge heights')
    pause(1)
end
return
%% Edges
fname = sprintf('edges_%s',datestr(now,'ddmmyy_hhMMss'))
save(fname,'edges')

