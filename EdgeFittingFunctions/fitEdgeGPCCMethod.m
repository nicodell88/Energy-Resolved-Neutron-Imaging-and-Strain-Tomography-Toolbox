function [deltaD,sigma,TrFit,fitinfo,opts] = fitEdgeGPCCMethod(Tr1,Tr2,tof,opts)
%fitEdgeGPCCMethod fits a bragg-edge using the method presented in: 
% TODO: add reference to ArXiv etc
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
%       opts.ns     :   Number of samples to use in MC step.
%       opts.nx      :   Number test points.
% Outputs:
%   - edgePos is the location of the braggEdge
%   - sigma is the estimated standard deviation
%   - TrFit is is the Bragg edge model evaluated at tof
%   - fitinfo contains additional information about the quality of the fit
%       fitinfo.lengthscale     : the lengthscale used
%
%See also fitEdges.

%
% Copyright (C) 2020 The University of Newcastle, Australia
% Authors:
%   Nicholas O'Dell <Nicholas.Odell@newcastle.edu.au>
% Last modified: 23/06/2020
% This program is licensed under GNU GPLv3, see LICENSE for more details.

sig_f   = 1;
l       = 1e-4;     % GP lengthscale
ns      = 3000;     % number of samples to draw to compute variance
nx      = 2500;     % number of points to predict at
covfunc = 'm52';     % default covariance function
m_basis = 1000;     % default number of basis functions for hp optimisation

%GP
if isfield(opts,'covfunc')
    covfunc = opts.covfunc;
end
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



if isfield(opts,'optimiseHP')
    optimiseHP = opts.optimiseHP;
end
if isfield(opts,'save_samples')
    save_samples = opts.save_samples;
else
    save_samples = false;
end


Tr1 = Tr1 - mean(Tr1(1:20));
Tr2 = Tr2 - mean(Tr2(1:20));
M = mean(Tr2(end-20:end))/mean(Tr1(end-20:end));
Tr1 = Tr1*M;

Tr1 = Tr1(:);
Tr2 = Tr2(:);
tof = tof(:);
%% use GP to smooth and get derivative of both Tr1 and Tr2

y1 = Tr1;
y2 = Tr2;
x       = tof(:);
xt = linspace(tof(opts.startIdx(2)),tof(opts.endIdx(1)),nx)';
%% identify sigm for both d0 and d
sig_m1 = mean(std(Tr1(end-40:end)),std(Tr1(1:40)));
sig_m2 = std(Tr2(end-40:end));
%%

% sig_m = std([Tr(opts.endIdx(1):opts.endIdx(2)) - g2(tof(opts.endIdx(1):opts.endIdx(2))),...
%     Tr(opts.startIdx(1):opts.startIdx(2))-g1(tof(opts.startIdx(1):opts.startIdx(2)))]);

% dt = tof(opts.fitIdx(2))-tof(opts.fitIdx(2)-1);

% sig_m = 2e-2;


optimiseHP = 'none';  % if true optimises the lengthscale for each Transmission spectra


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
Phi = Phi;
[~,m]= size(Phi);
Gamma1 = [Phi./sig_m1;diag(1./sqrt(SLambda))];
Gamma2 = [Phi./sig_m2;diag(1./sqrt(SLambda))];
R1 = triu(qr(Gamma1));
R2 = triu(qr(Gamma2));
CZ1 = R1(1:m,1:m);
CZ2 = R2(1:m,1:m);
LSoptsT.TRANSA = true; LSoptsT.UT = true; LSopts.TRANSA = false; LSopts.UT = true;
v_basis1 = (linsolve(CZ1,linsolve(CZ1,Phi'*(y1./sig_m1.^2),LSoptsT),LSopts));
v_basis2 = (linsolve(CZ2,linsolve(CZ2,Phi'*(y2./sig_m2.^2),LSoptsT),LSopts));

festp1 = Phi_p*v_basis1;  %
g1 = dPhi_T * v_basis1;

festp2 = Phi_p*v_basis2;  %
dTr2 = dPhi_T * v_basis2;

% alternatively
gamma = linsolve(CZ1,eye(m),LSoptsT);
SigV = triu(qr(gamma));
sv = v_basis1  + SigV.'*randn(m,ns);      % sample the coefficients

sv2 = v_basis2 + SigV.'*randn(m,ns);      % sample the coefficients
%% sample the derivative of the projection and use the mean of the d0
sg  = dPhi_T * sv;
sg2 = dPhi_T * sv2;
%% calculate normalised cross correlation for all samples
XC = nan(nx*2-1,ns);
for i =1:ns
    [XC(:,i),lags] = xcorr(sg(:,i),dTr2,[],'normalized');
    %     [XC(:,i),lags] = xcorr(sg(:,i),sg2(:,i),[],'normalized');
end

dt = xt(2)-xt(1);
%% take max
[~,I] = max(XC,[],1);
%% calculate sample mean and standard deviation
d_all = lags(I) * dt;

deltaD = mean(d_all);
sigma  = std(d_all);
TrFit  = festp1;

fitinfo.lengthscale = l;
end