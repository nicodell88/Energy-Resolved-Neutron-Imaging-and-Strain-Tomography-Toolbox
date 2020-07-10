%
% runGetD0.m
%
% Copyright (C) 2020 The University of Newcastle, Australia
% Authors:
%   Nicholas O'Dell <Nicholas.Odell@newcastle.edu.au>
% Last modified: 10/07/2020
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
OBfile = '/Users/megatron/Dropbox/Strain Tomography 2020/J-PARC 2019 Working Folder/data/open_beam/open_beam_11ksec.mat';
D0file = '/Users/megatron/Dropbox/Strain Tomography 2020/J-PARC 2019 Working Folder/data/wedge/wedge_11ksec.mat';

[Tr,tof] = getD0(OBfile,D0file);
%% Save results
%  save('d0Profile','Tr')
%% 
figure(2)
clf
plot(tof,Tr)
title('D0 profile')
ylabel('Transmission Intensity')
xlabel('Time-of-flight - [seconds]')