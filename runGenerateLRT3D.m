% runGenerateLRT3D
%
% Copyright (C) 2020 The University of Newcastle, Australia
% Authors:
%   Nicholas O'Dell <Nicholas.Odell@newcastle.edu.au>
% Last modified: 15/07/2020
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
ProjIdx = 0:69;  %Projections to process
 ProjFileSpec    = ...'/Users/megatron/Dropbox/Strain Tomography 2020/J-PARC 2020 Working Folder/Data/OpenBeam/OpenBeamImg__proj001_idx.mat';
                      '/Users/megatron/Dropbox/Strain Tomography 2020/J-PARC 2019 Working Folder/data/cube_plug/Cube1__proj%d_idx.mat';
OBfile          = ...'/Users/megatron/Dropbox/J-Parc2018 Working Folder/data/ProcessedWithToolbox/preprocessed/OpenBeam_proj_003.mat';
                     '/Users/megatron/Dropbox/Strain Tomography 2020/J-PARC 2019 Working Folder/data/open_beam/open_beam_11ksec.mat';       

STL_file         = './Alignment/cubePlugCorrected.stl';
AlignResult      = load('./SampleAlignmentResults.mat');

nanThresh = 0.8;
%% Tuning parametes for averaging
nPix = 32;

[Tr,wl,tof,entry,exit,L,nhat,yInds,nSegs] = GenerateLRT3D(ProjIdx,ProjFileSpec,OBfile,STL_file,'mm',AlignResult.mask,AlignResult.rODn,AlignResult.Rno,nPix,nanThresh);
%                                                        (ProjIdx,ProjFileSpec,OBfile,stlFile ,stlUnit,mask_all,rODn_all,Rno_all,nPix,nanThresh,varargin)

