% runProcessLRT2Dexample
%
% Copyright (C) 2020 The University of Newcastle, Australia
% Authors:
%   Nicholas O'Dell <Nicholas.Odell@newcastle.edu.au>
% Last modified: 06/07/2020
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
delete(findall(0,'tag','TMWWaitbar'));
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

%% Inputs
ProjIdx = 0:49;  %Projections to process
ProjFileSpec    = '/Users/megatron/Dropbox/J-Parc2018 Working Folder/data/ProcessedWithToolbox/preprocessed/Sample2018_proj_%03d.mat';
OBfile          = '/Users/megatron/Dropbox/J-Parc2018 Working Folder/data/ProcessedWithToolbox/preprocessed/OpenBeam_proj_003.mat';

sample      = load('/Users/megatron/Dropbox/J-Parc2018 Working Folder/data/ProcessedWithToolbox/BoundaryRP.mat');
AlignResult = load('/Users/megatron/Dropbox/J-Parc2018 Working Folder/data/ProcessedWithToolbox/RPmask.mat');

%% Tuning parametes for averaging
colRange   = [20,250];
nRow       = 4;

[Tr,wl,tof,entry,exit,L,nhat,yInds,nSegs] = GenerateLRT2D(ProjIdx,ProjFileSpec,OBfile,sample.BoundaryRP,AlignResult.mask,AlignResult.rODn,AlignResult.Rno,colRange,nRow);


