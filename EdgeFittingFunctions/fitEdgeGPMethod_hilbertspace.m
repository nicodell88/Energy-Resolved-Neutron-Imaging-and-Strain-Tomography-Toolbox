function [edgePos,sigma,TrFit,fitinfo] = fitEdgeGPMethod_hilbertspace(Tr,tof,opts)
%fitEdgeGPMethod fits a bragg-edge using the method presented in:
%   TODO: Insert Bib entry
%   TODO: Insert arxiv link when paper is written.
%
% Inputs:
%   - Tr is a 1xN double containing the normalised transmisssion curve
%   for a single projection
%   - tof is an 1xN array of wave-lengths or time-of-flight.
%   - options is a structure containing
%       opts.a00    :   Initial guess
%       opts.b00    :   Initial guess
%       opts.a_hkl0 :   Initial guess
%       opts.b_hkl0 :   Initial guess
%       opts.sig_f  :   Squared-Exponential Kernel Hyperparameter, output variance
%       opts.l      :   Squared-Exponential Kernel Hyperparameter, lengthscale
%
% Outputs:
%   - edgePos is the location of the braggEdge
%   - sigma is the estimated standard deviation
%   - TrFit is is the Bragg edge model evaluated at tof
%   - fitinfo contains additional information about the quality of the fit
%       fitinfo.lengthscale     : the lengthscale used
%       fitinfo.std_residual    : the standard deviaton of the residual
%       fitinfo.rms_residual    : the root mean square of the residual
%       fitinfo.fitqual         : an estimate of the fit quality given by
%                                   the radio sig_m / std(residual)
%       fitinfo.widthathalfheight: the peak width at half height
%       fitinfo.edgeshape       : the estimated edge shape
%       fitinfo.gradedgeshape   : the gradient of the estimated edge shape
%       fitinfo.x_grad          : the x (tof/wavelength) that the gradient
%                                   is calculated at
%
%
% Copyright (C) 2020 The University of Newcastle, Australia
% Authors:
%   Johannes Hendriks <Johannes.Hendriks@newcastle.edu.au>
% Last modified: 25/03/2020
% This program is licensed under GNU GPLv3, see LICENSE for more details.

%% least squares fitting options
optionsFit              = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');
optionsFit.Algorithm    = 'Levenberg-Marquardt';
optionsFit.Jacobian     = 'off';
optionsFit.Display      = 'off';
%% Initial guess
a00     = 0.5;
b00     = 0.5;
a_hkl0  = 0.5;
b_hkl0  = 0.5;
sig_f   = 1;    
l       = 1e-4;     % GP lengthscale
ns      = 3000;     % number of samples to draw to compute variance
nx      = 2500;     % number of points to predict at
optimiseHP = false;  % if true optimises the lengthscale for each Transmission spectra
covfunc = 'se';     % default covariance function
m_basis = 1000;     % default number of basis functions for hp optimisation

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


%GP
if isfield(opts,'sig_f')
    sig_f = opts.sig_f;
end
if isfield(opts,'l')        % default lengthscales
    l = opts.l;
else
    if strcmpi(covfunc,'se') 
        l = 0.03;
    end
    if strcmpi(covfunc,'m32') 
        l = 0.18;
    end
    if strcmpi(covfunc,'m52')
        l = 0.09;
    end
end
if isfield(opts,'numberbasis')
    m_basis = opts.numberbasis;
end
if isfield(opts,'ns')
    ns = opts.ns;
end
if isfield(opts,'nx')
    nx = opts.nx;
end
if isfield(opts,'GPscheme')
   GPscheme = opts.GPscheme; 
   if ~strcmpi(GPscheme,'hilbertspace')
       error('Should not have gotten here')
    end
end

if isfield(opts,'covfunc')
    covfunc = opts.covfunc;
end

if isfield(opts,'optimiseHP')
    optimiseHP = opts.optimiseHP;
end
if isfield(opts,'save_samples')
    save_samples = opts.save_samples;
else
    save_samples = false;
end


%% Fit edge
%% 1) fit to the far right of the edge where B = 1, so only fit exp([-(a0+b0.*t)])
fit1 = @(p,x) exp(-(p(1) + p(2).*x));
[p,~,~,~,~,~,~] = lsqcurvefit(fit1,[a00;b00],tof(opts.endIdx(1):opts.endIdx(2)),Tr(opts.endIdx(1):opts.endIdx(2)),[],[],optionsFit);
a0 = p(1); b0 = p(2);
%% 2) fit to the far left of the edge where B = 0;
fit2 = @(p,x) exp(-(a0 + b0.*x)).*exp(-(p(1)+p(2).*x));
[p,~,~,~,~,~,~] = lsqcurvefit(fit2,[a_hkl0;b_hkl0],tof(opts.startIdx(1):opts.startIdx(2)),Tr(opts.startIdx(1):opts.startIdx(2)),[],[],optionsFit);
a_hkl = p(1); b_hkl = p(2);
%% 3) fit the transition function as a GP
g1 = @(x) exp(-(a0 + b0.*x)).*exp(-(a_hkl+b_hkl.*x));
g2 = @(x) exp(-(a0 + b0.*x));

y = (Tr - g1(tof)).';
x = tof.';
sig_m = std([Tr(opts.endIdx(1):opts.endIdx(2)) - g2(tof(opts.endIdx(1):opts.endIdx(2))),...
    Tr(opts.startIdx(1):opts.startIdx(2))-g1(tof(opts.startIdx(1):opts.startIdx(2)))]);
% nh = nx;
% nx = length(tof); 
xt = linspace(tof(opts.startIdx(2)),tof(opts.endIdx(1)),nx)';
% xt_interp = linspace(tof(opts.startIdx(2)),tof(opts.endIdx(1)),nh)';


% Hyperparameters
if optimiseHP
    if strcmpi(covfunc,'se')
        fminopts = optimoptions('fminunc','SpecifyObjectiveGradient',true,'display','none');
        nlM = @(l) LogMarginalSE(l,x,y,g1(x),g2(x),sig_m);
        logl = fminunc(nlM,0,fminopts);
        l = max((tof(2)-tof(1))*60,exp(logl));      % ensure a sensible result
    elseif strcmpi(covfunc,'m32')
        fminopts = optimoptions('fminunc','SpecifyObjectiveGradient',true,'display','none');
        nlM = @(l) LogMarginalM32(l,x,y,g1(x),g2(x),sig_m);
        logl = fminunc(nlM,0,fminopts);
        l = max((tof(2)-tof(1))*10,exp(logl));      % ensure a sensible result
    elseif strcmpi(covfunc,'m52')
        fminopts = optimoptions('fminunc','SpecifyObjectiveGradient',false,'display','none');
        nlM = @(l) LogMarginalM52(l,x,y,g1(x),g2(x),sig_m);
        logl = fminunc(nlM,0,fminopts);
        l = max((tof(2)-tof(1))*30,exp(logl));      % ensure a sensible result
    end
end


if strcmpi(covfunc,'se')
    dlambda = 3.5/l/m_basis;    % 3.5 sigma coverage
    L = max(pi/2/dlambda,tof(end)-tof(1));
    [Phi,~,SLambda,~, dPhi_T] = hilbert_approxSE(l,sig_f,m_basis,L,xt,x);
elseif strcmpi(covfunc,'m32')
    dlambda = 30/l/m_basis;    
    L = max(pi/2/dlambda,tof(end)-tof(1));
    [Phi,~,SLambda,~, dPhi_T] = hilbert_approxM32(l,sig_f,m_basis,L,xt,x);  

elseif strcmpi(covfunc,'m52')
    dlambda = 11/l/m_basis;    
    L = max(pi/2/dlambda,tof(end)-tof(1));
    [Phi,~,SLambda,~, dPhi_T] = hilbert_approxM52(l,sig_f,m_basis,L,xt,x);  

else 
    error('Invalid covariance function.')
end

Phi_p = Phi;
Phi = Phi .* (g2(x) - g1(x));
[~,m]= size(Phi);
Gamma = [Phi./sig_m;diag(1./sqrt(SLambda))];
R = triu(qr(Gamma));
CZ = R(1:m,1:m);
LSoptsT.TRANSA = true; LSoptsT.UT = true; LSopts.TRANSA = false; LSopts.UT = true;
v_basis = (linsolve(CZ,linsolve(CZ,Phi'*(y./sig_m.^2),LSoptsT),LSopts));
festp = Phi_p*v_basis;  %
g = dPhi_T * v_basis;

% alternatively
gamma = linsolve(CZ,eye(m),LSoptsT);
SigV = triu(qr(gamma));
sv = v_basis + SigV.'*randn(m,ns);      % sample the coefficients
sg = dPhi_T * sv;


[~,Is] = max([g sg]);
sLams = xt(Is);


%% Collect Results
edgePos = mean(sLams);
sigma = std(sLams);
TrFit = exp(-a0-b0*tof).*...
	(exp(-a_hkl-b_hkl*tof) + (1-exp(-a_hkl - b_hkl*tof)) .*festp.');

% check fit quality
fitqual = sig_m/std(Tr-TrFit);
if fitqual > 2
    warning('The ratio of sig_m/std(residual) is high ( %s), indicating that the data may have been overfit. Consider increasing the lengthscale',num2str(fitqual))
end
if fitqual < 0.1
    warning('The ratio of sig_m/std(residual) is low ( %s), indicating that the data may have been overfit. Consider increasing the lengthscale',num2str(fitqual))
end


fitinfo.lengthscale = l;                            % store the lengthscale used
fitinfo.std_residual = std(Tr-TrFit);               % standard deviation of the residual
fitinfo.rms_residual = sqrt(mean((Tr-TrFit).^2));   % root mean square of hte residual
fitinfo.fitqual = fitqual;
fitinfo.edgeshape = festp;
fitinfo.gradedgeshape = g;
fitinfo.x_grad = xt;

if save_samples
    fitinfo.lambda_samples = sLams;    
end

half_height = max(g)/2;

[xi,~] = polyxpoly(xt,g,[xt(1);xt(end)],[half_height;half_height]);
if length(xi) == 2
    widthathalfheight = max(xi) - min(xi);
else
    widthathalfheight = nan;
end
fitinfo.widthathalfheight = widthathalfheight;





end