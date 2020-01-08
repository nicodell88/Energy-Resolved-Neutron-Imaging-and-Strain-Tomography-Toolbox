function [d_cell,std_cell,TrFit_cell] = fitEdges(Tr,tof,opts)
%FITEDGES Fits Bragg edges to input data.
%   [d_cell,std_cell,TrFit_cell] = fitEdges(Tr,tof,opts)
%   Inputs:
%       - Tr is a cell-array of normalised transmission intensity data. Each cell
%       is a projection.
%       - tof is an array of wave-lengths of time-of-flight, whos length matches
%       the width of each array in Tr.
%       - options is a structure containing
%           opts.method      :  a char array defining the fitting method
%                               ('Santisteban2001', 'Tremsin2011', 'Hendriks2020')
%           opts.startRange  :  a 2 element vector containing the start
%                               and end range for fitting the left side of the
%                               Bragg-Edge [start end] with the same units as
%                               tof
%           opts.endRange    :  a 2 element vector containing the start
%                               and end range for fitting the right side of
%                               the Bragg-Edge [start end] with the same
%                               units as tof
%           opts.plot        :  a logical flag. When true this will plot
%                               each Bragg edge as it is fitted.
%   Outputs:
%       - d_cell is a cell array containing the bragg edge location for
%       each projection.
%       - std_cell is a cell array containing corresponding standard
%       deviation estimates for each result in d_cell.
%       - TrFit is is the Bragg edge model for each fit evaluated at tof
%
% Copyright (C) 2020 The University of Newcastle, Australia
% Authors:
%   Nicholas O'Dell <Nicholas.Odell@newcastle.edu.au>
% Last modified: 08/01/2020
% This program is licensed under GNU GPLv3, see LICENSE for more details.

%% TODO: alow opts to contain initial guess for curve elements

%% Inspect Data
assert(numel(tof)==length(tof),'Expected tof to be a Nx1 vector');
%TODO: check that tr is a cell array
%TODO: test that length(tof) matches dimension 2 of each cell in Tr

tof = tof(:).';
%% Fill opts structure
if nargin==3
    %% Method
    if isfield(opts,'method')
        switch lower(opts.method)
            case 'santisteban2001'
                edgeFit = @(tr,wl,op)fitEdgeSantisteban2001(tr,wl,op);
            case 'tremsin2011'
                edgeFit = @(tr,wl,op)fitEdgeTremsin2011(tr,wl,op);
            case 'hendriks2020'
                edgeFit = @(tr,wl,op)fitEdgeSantisteban2001(tr,wl,op);
            otherwise
                error('%s is not a valid edge fitting method, see help fitEdges',opts.method);
        end
    else
        opts.method = 'Santisteban2001';
        edgeFit = @(tr,wl,op)fitEdgeSantisteban2001(tr,wl,op);
    end
    %% left and right ranges
    if isfield(opts,'startRange')
        assert(numel(opts.startRange)==2,'Expected opts.startRange to contain two elements.')
        [~,opts.startIdx] = min((tof(:).'-opts.startRange(:)).^2,[],2);
    else
        ntof = length(tof);
        opts.startIdx    = [1 round(0.25*ntof)];
    end
    
    if isfield(opts,'endRange')
        assert(numel(opts.endRange)==2,'Expected opts.endRange to contain two elements.')
        [~,opts.endIdx] = min((tof(:).'-opts.endRange(:)).^2,[],2);
    else
        ntof = length(tof);
        opts.endIdx      = [round(0.75*ntof) ntof];
    end
else
    opts.method = 'Santisteban2001';
    edgeFit = @(tr,wl,op)fitEdgeSantisteban2001(tr,wl,op);
    ntof = length(tof);
    opts.startIdx    = [1 round(0.25*ntof)];
    opts.endIdx      = [round(0.75*ntof) ntof];
    opts.plot   = false;
end

%% fitEdge
np = numel(Tr); %number of projections

% Initialise cells for the fitted d-spacings and standard deviations
d_cell = cellfun(@(x) NaN(size(x,1),1), Tr, 'UniformOutput',false);
std_cell = cellfun(@(x) NaN(size(x,1),1), Tr, 'UniformOutput',false);
TrFit_cell = cellfun(@(x) NaN(size(x)), Tr, 'UniformOutput',false);

if opts.plot
   figure
   Hdata    = plot(tof,nan(size(tof)),'.');
   hold on
   Hfit     = plot(tof,nan(size(tof)),'.');
   Htitle   = title('Projection');
   grid minor
   ylim([0 1])
end

wh = updateWaitbar();
% Loop over projections
for k = 1:np
  [wh,flag] = updateWaitbar(k,np,wh);
  if flag
      break
  end
   %Loop over measurements in this projection
   for i = 1:size(Tr{k},1)
       % Call edge fitting function
       [d_cell{k}(i),std_cell{k}(i),TrFit_cell{k}(i,:)] = edgeFit(Tr{k}(i,:),tof,opts);
       % Plot Results
       if opts.plot
           msg = sprintf('Projection %d, Measurement %d',k,i);
           Htitle.String = msg;
           Hdata.YData = Tr{k}(i,:);
           Hfit.YData  = TrFit_cell{k}(i,:);
           drawnow
       end
   end
end
delete(wh)
end

