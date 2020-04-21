function [d_cell,std_cell,TrFit_cell,fitinfo_cell,opts] = fitEdges(Tr,tof,opts,d0)
%FITEDGES Fits Bragg edges to input data.
%   [d_cell,std_cell,TrFit_cell] = fitEdges(Tr,tof,opts)
%   Inputs:
%       - Tr is a cell-array of normalised transmission intensity data. Each cell
%       is a projection.
%       - tof is an array of wave-lengths of time-of-flight, whos length matches
%       the width of each array in Tr.
%       - options is a structure containing
%           opts.method      :  a char array defining the fitting method
%                               ('attenuation', '5param', 'gp')
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
%           opts. ...        :  Other options specific to fitEdge...
%   Outputs:
%       - d_cell is a cell array containing the bragg edge location for
%       each projection.
%       - std_cell is a cell array containing corresponding standard
%       deviation estimates for each result in d_cell.
%       - TrFit is is the Bragg edge model for each fit evaluated at tof
%       - fitinfo_cell contains additional fit information
%
% Copyright (C) 2020 The University of Newcastle, Australia
% Authors:
%   Nicholas O'Dell <Nicholas.Odell@newcastle.edu.au>
%   Johannes Hendriks <Johannes.Hendriks@newcastle.edu.au>
% Last modified: 21/04/2020
% This program is licensed under GNU GPLv3, see LICENSE for more details.
TBdir = fileparts(mfilename('fullpath'));
addpath(fullfile(TBdir,'EdgeFittingFunctions'));
addpath(genpath(fullfile(TBdir,'utility_functions')))
%% Inspect Data
assert(numel(tof)==length(tof),'Expected tof to be a Nx1 vector');
assert(iscell(Tr),'Expected ''Tr'' to be a cell array');
[~,Trncols] = cellfun(@size, Tr);
assert(all(Trncols == length(tof)),'Expected the number of columns in each cell of Tr to be equal to the length of tof');

tof = tof(:).';
%% Fill opts structure
if nargin>=3
    %% Method
    if isfield(opts,'method')
        switch lower(opts.method)
            case 'attenuation'
                edgeFit = @(tr,wl,op)fitEdgeAttenuationMethod(tr,wl,op);
            case 'attenuationalt'
                edgeFit = @(tr,wl,op)fitEdgeAttenuationMethodAlt(tr,wl,op);
            case '5param'
                edgeFit = @(tr,wl,op)fitEdge5ParamMethod(tr,wl,op);
            case 'gp'
                if isfield(opts,'GPscheme')
                    if ~isfield(opts,'covfunc')
                        opts.covfunc = 'SE';
                    end
                    if strcmpi(opts.GPscheme,'hilbertspace') || ~strcmpi(opts.covfunc,'SE')
                        if ~strcmpi(opts.covfunc,'SE') && ~strcmpi(opts.GPscheme,'hilbertspace')
                            warning('%s covariance func can only be used with hilbertspace scheme.',opts.GPscheme)
                        end
                        edgeFit = @(tr,wl,op)fitEdgeGPMethod_hilbertspace(tr,wl,op);
                    elseif strcmpi(opts.GPscheme,'full') || strcmpi(opts.GPscheme,'interp')
                        edgeFit = @(tr,wl,op)fitEdgeGPMethod(tr,wl,op);
                    else
                        error('%s is not a valid GP scheme, see help fitEdges',opts.GPscheme);
                    end
                else
                    edgeFit = @(tr,wl,op)fitEdgeGPMethod(tr,wl,op);
                end
            case 'crosscorr'
                edgeFit = @(tr,wl,op)crossCorrMethod(tr,d0,wl,op);
            otherwise
                error('%s is not a valid edge fitting method, see help fitEdges',opts.method);
        end
    else
        opts.method = 'GP';
        opts.GPscheme = 'hilbertspace';
        opts.covfunc = 'M52';
        edgeFit = @(tr,wl,op)fitEdgeGPMethod_hilbertspace(tr,wl,op);
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
    
    if isfield(opts,'range')
        assert(numel(opts.range)==2,'Expected opts.range to contain two elements.')
        [~,opts.rangeIdx] = min((tof(:).'-opts.range(:)).^2,[],2);
    else
        ntof = length(tof);
        opts.rangeIdx     = [round(0.35*ntof) round(0.65*ntof)];
    end
    if ~isfield(opts,'plot')
        opts.plot =false;
    end
    if ~isfield(opts,'Par')
        opts.Par = false;
    end
else
    opts.method = 'attenuation';
    edgeFit = @(tr,wl,op)fitEdgeAttenuation(tr,wl,op);
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
[~,~,~,fitInfoTemplate] = edgeFit(Tr{1}(1,:),tof,opts);
fitinfo_cell= cellfun(@(x) repmat(fitInfoTemplate,size(x,1),1), Tr, 'UniformOutput',false);

if opts.plot && ~opts.Par
    if strcmpi(opts.method,'5param')
        plot_idx = opts.rangeIdx(1):opts.rangeIdx(2);
    elseif strcmpi(opts.method,'crosscorr')
    else
        plot_idx = opts.startIdx(1):opts.endIdx(end);
    end
    Hfig = figure;
    Hdata    = plot(tof(plot_idx),nan(size(tof(plot_idx))),'.');
    hold on
    Hfit     = plot(tof(plot_idx),nan(size(tof(plot_idx))),'--');
    Htitle   = title('Projection','Interpreter','Latex');
    xlabel('[Wave Length] or \{Time of Flight\} - [\AA] or \{s\}','Interpreter','Latex')
    legend('Tr','Edge fit')
    grid minor
    ylim([0 1])
    hold off
end


[Cnrows,~] = cellfun(@size, Tr);
nAll = sum(Cnrows);
%% fit Edges
if ~opts.Par %no parfor
    % Initialise Waitbar
    wh = updateWaitbar();
    iter = 0;
    % Loop over projections
    for k = 1:np
        %Loop over measurements in this projection
        for i = 1:size(Tr{k},1)
            iter = iter+1;
            [wh,flag] = updateWaitbar(iter,nAll,wh);
            if flag
                return
            end
            % Call edge fitting function
            try
                [d_cell{k}(i),std_cell{k}(i),TrFit_cell{k}(i,:),fitinfo_cell{k}(i)] = edgeFit(Tr{k}(i,:),tof,opts);
            catch e
                delete(wh)
                if opts.plot
                    close(Hfig);
                end
%                 error(e.message);
                 rethrow(e)
            end
            % Plot Results
            if opts.plot
                msg = sprintf('Projection %d, Measurement %d',k,i);
                Htitle.String = msg;
                Hdata.YData = Tr{k}(i,plot_idx);
                Hfit.YData  = TrFit_cell{k}(i,plot_idx);
                drawnow
            end
        end
    end
    delete(wh)
    if opts.plot
        close(Hfig);
    end
elseif opts.Par &&  np==1 %par for single proj
        k=1;%only 1 projection
        %Loop over measurements in this projection
        nEdge = size(Tr{k},1);
        
        dTemp = nan(nEdge,1);
        stdTemp = nan(nEdge,1);
        trFitTemp = nan(size(Tr{k}));
        fitInfoTemp = repmat(fitInfoTemplate,nEdge,1);
        
        % Setup par for waitbar
        ppm = ParforProgressbar(nAll,'title','Fitting Bragg Edges'); 
        
        parfor i = 1:nEdge
            ppm.increment();
            % Call edge fitting function
            try
%                     [dTemp(i),stdTemp(i),trFitTemp(i,:),~] = fitEdgeAttenuationMethod(Tr{k}(i,:),tof,opts);
                  [dTemp(i),stdTemp(i),trFitTemp(i,:),fitInfoTemp(i)] = edgeFit(Tr{k}(i,:),tof,opts);
%                 fitInfoTemp(i) = fitInfoTemptmp;
            catch e
%                 error(e.message);
            delete(ppm);
            rethrow(e)
            end
            
        end
        delete(ppm);
        d_cell{k}       = dTemp;
        std_cell{k}     = stdTemp;
        TrFit_cell{k}   = trFitTemp;
        fitinfo_cell{k} = fitInfoTemp;
else %par for multipl proj
    ppm = ParforProgressbar(nAll,'title','Fitting Bragg Edges'); 
    parfor k = 1:np
        %Loop over measurements in this projection
        nEdge = size(Tr{k},1);
        
        dTemp = nan(nEdge,1);
        stdTemp = nan(nEdge,1);
        trFitTemp = nan(size(Tr{k}));
        fitInfoTemp = repmat(fitInfoTemplate,nEdge,1);
        for i = 1:nEdge
            ppm.increment();
            % Call edge fitting function
            try
%                 [dTemp(i),stdTemp(i),trFitTemp(i,:),~] = fitEdgeAttenuationMethod(Tr{k}(i,:),tof,opts);
                [dTemp(i),stdTemp(i),trFitTemp(i,:),fitInfoTemp(i)] = edgeFit(Tr{k}(i,:),tof,opts);
%                 fitInfoTemp(i) = fitInfoTemptmp;
            catch e
%                 error(e.message);
            delete(ppm);
            rethrow(e)
            end
            
        end
        d_cell{k}       = dTemp;
        std_cell{k}     = stdTemp;
        TrFit_cell{k}   = trFitTemp;
        fitinfo_cell{k} = fitInfoTemp;
    end
    delete(ppm);
end%if

end%function

