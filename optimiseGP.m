function [lengthscale,opts] = optimiseGP(Tr,tof,opts)
%FITEDGES Optimise the lengthscale used for Gaussian process regression of
%the bragg edges
%   [lengthscale,opts] = optimiseGP(Tr,tof,opts)
%   Inputs:
%       - Tr is a cell-array of normalised transmission intensity data. Each cell
%       is a projection.
%       - tof is an array of wave-lengths of time-of-flight, whos length matches
%       the width of each array in Tr.
%       - options is a structure containing
%           opts.method      :  must be set to GP
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
%       - lengthscale: the optimum lengthscale for the particular scheme
%       and covariance function specified in opts
%       - opts: the options structure with the lengthscale set and
%       optimiseHP set to false so that optimisation wont be rerun
%
% Copyright (C) 2020 The University of Newcastle, Australia
% Authors:
%   Johannes Hendriks <Johannes.Hendriks@newcastle.edu.au>
% Last modified: 02/04/2020
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
if nargin==3
    %% Method
    if isfield(opts,'method')
        if ~any(strcmpi(opts.method,{'gp','gpcc2'}))
            error('Can only optimise HP if using GP method')
        else
%             edgeFit = @(tr,wl,op)fitEdgeGPMethod(tr,wl,op);
        end
    else
        warning('Method not specified, setting method to GP.')
%         edgeFit = @(tr,wl,op)fitEdgeGPMethod(tr,wl,op);
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
else
    warning('Method not specified, Setting method to GP')
    opts.method = 'gp';
%     edgeFit = @(tr,wl,op)fitEdgeGPMethod(tr,wl,op);
    ntof = length(tof);
    opts.startIdx    = [1 round(0.25*ntof)];
    opts.endIdx      = [round(0.75*ntof) ntof];
    opts.plot   = false;
end

if isfield(opts,'covfunc')
    if strcmpi(opts.covfunc,'se')
        logpdffunc = @(logl,X,Y,g1,g2,sig_e) LogMarginalSE(logl,X,Y,g1,g2,sig_e);
    elseif strcmpi(opts.covfunc,'m32')
        logpdffunc = @(logl,X,Y,g1,g2,sig_e) LogMarginalM32(logl,X,Y,g1,g2,sig_e);
    elseif strcmpi(opts.covfunc,'m52')
        logpdffunc = @(logl,X,Y,g1,g2,sig_e) LogMarginalM52(logl,X,Y,g1,g2,sig_e);
    else
        error('Invalid covariance function specified') 
    end
else
    warning('Covfunc not specified, setting covariance function to squared-exponential.')
    logpdffunc = @(logl,X,Y,g1,g2,sig_e) LogMarginalSE(logl,X,Y,g1,g2,sig_e);
end

%% package data for optimisatoin
np = numel(Tr); %number of projections
% least squares fitting options
optionsFit              = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');
optionsFit.Algorithm    = 'Levenberg-Marquardt';
optionsFit.Jacobian     = 'off';
optionsFit.Display      = 'off';
a00     = 0.5;
b00     = 0.5;
a_hkl0  = 0.5;
b_hkl0  = 0.5;
if isfield(opts,'a00')
    a00 = opts.a00;
end
if isfield(opts,'b00')
    b00 = opts.b00;
end
if isfield(opts,'a_hkl0')
    a_hkl0 = opts.a_hkl0;
end
if isfield(opts,'b_hkl0')
    b_hkl0 = opts.b_hkl0;
end

tr = [];
g1 = [];
g2 = [];
Y = [];
Sig_m = [];
% Loop over projections
for k = 1:np
    %Loop over measurements in this projection
    for i = 1:size(Tr{k},1)
        try
        fit1 = @(p,x) exp(-(p(1) + p(2).*x));
        [p,~,~,~,~,~,~] = lsqcurvefit(fit1,[a00;b00],tof(opts.endIdx(1):opts.endIdx(2)),Tr{k}(i,opts.endIdx(1):opts.endIdx(2)),[],[],optionsFit);
        a0 = p(1); b0 = p(2);
        %% 2) fit to the far left of the edge where B = 0;
        fit2 = @(p,x) exp(-(a0 + b0.*x)).*exp(-(p(1)+p(2).*x));
        [p,~,~,~,~,~,~] = lsqcurvefit(fit2,[a_hkl0;b_hkl0],tof(opts.startIdx(1):opts.startIdx(2)),Tr{k}(i,opts.startIdx(1):opts.startIdx(2)),[],[],optionsFit);
        a_hkl = p(1); b_hkl = p(2);

        g1 = [g1,exp(-(a0 + b0.*tof.')).*exp(-(a_hkl+b_hkl.*tof.'))];
        g2 = [g2,exp(-(a0 + b0.*tof.'))];
        Y = [Y,Tr{k}(i,:).' - g1(:,end)];
        tr = [tr;Tr{k}(i,:)];
        
        sig_m = std([Tr{k}(i,opts.endIdx(1):opts.endIdx(2)).' - g2(opts.endIdx(1):opts.endIdx(2),end);...
                Tr{k}(i,opts.startIdx(1):opts.startIdx(2)).'-g1(opts.startIdx(1):opts.startIdx(2),end)]);
        Sig_m = [Sig_m;sig_m];
        catch
           warning(['Was unable to use data from projection ', k, ' edge ', i, ' for optimisatoin due to failure to fit the exponentials']) 
        end
    end
end

if isempty(Y)
    warning('Failed to run optimisation')
    lengthscale = NaN;
    return
end

[~,c] = size(Y);
if c > 200
    warning(['About to run optimisation on data from ', num2str(c), ' Bragg-edges.',...
 'This could be very slow, consider passing in a subset of the data.'])
end

disp('Optimising lengthscale')
if ~strcmpi(opts.covfunc,'m52')
    fminopts = optimoptions('fminunc','SpecifyObjectiveGradient',true,'display','iter');
    nlM = @(logl) sumlogpdf(logpdffunc,logl,tof.',Y,g1,g2,Sig_m);
else
    fminopts = optimoptions('fminunc','SpecifyObjectiveGradient',false,'display','iter');
    nlM = @(logl) sumlogpdfNG(logpdffunc,logl,tof.',Y,g1,g2,Sig_m);
end
logl = fminunc(nlM,0,fminopts);
lengthscale = max((tof(2)-tof(1))*10,exp(logl)); % ensure a sensible result


opts.l = lengthscale;

disp(['optimisation finished with l = ',num2str(lengthscale)]) 


end


function [logpdf, gradlogpdf] = sumlogpdf(logpdffunc,logl,X,Y,g1,g2,sig_e)
      [~,c] = size(Y);
      logpdf = 0;
      gradlogpdf=0;
      for i=1:c
          [ltmp,gtmp] = logpdffunc(logl,X,Y(:,i),g1(:,i),g2(:,i),sig_e(i));
          logpdf = logpdf + ltmp;
          gradlogpdf = gradlogpdf + gtmp;
          
      end
     
end

function [logpdf] = sumlogpdfNG(logpdffunc,logl,X,Y,g1,g2,sig_e)
      [~,c] = size(Y);
      logpdf = 0;
      for i=1:c
          [ltmp] = logpdffunc(logl,X,Y(:,i),g1(:,i),g2(:,i),sig_e(i));
          logpdf = logpdf + ltmp;
          
      end
     
end
