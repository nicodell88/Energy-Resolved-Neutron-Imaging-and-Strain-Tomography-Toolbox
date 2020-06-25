%
% runExampleGP.m
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
load CubeProj.mat
%% Plot Data
figure(1); clf;
plot(tof,Tr{3}(100,:))%arbitrary
xlabel('Time-Of-Flight - [seconds]')
ylabel('Normalised Transmission Intensity - [arbitrary units]')
grid minor

%% Subsample data for hyperparameter optimisation
opts = BraggOptions(tof,'gp','Par',true);
NBatch = 200;
[Cnrows,~] = cellfun(@size, Tr);

r = rand(NBatch,1);
idxP = randsample(numel(Tr),NBatch,'true');
idxE = round(r.*Cnrows(idxP));
extract = @(idxP_cell,idxE_cell) Tr{idxP_cell}(idxE_cell,:);
Tr_batch = {cell2mat(arrayfun(extract,idxP,idxE,'uniformoutput',false))};

[~,opts] = optimiseGP(Tr_batch,tof,opts);
%% FIT EDGES
[d_cell,std_cell,TrFit_cell,fitInfo_cell] = fitEdges(Tr,tof,opts);
