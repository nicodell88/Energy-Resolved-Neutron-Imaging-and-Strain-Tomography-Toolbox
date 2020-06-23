function [deltaD,sigma,TrFit,fitinfo,opts] = fitEdgeGPCC2(Tr1,Tr2,tof,opts)
%fitEdgeGPCC2 fits a bragg-edge using the method presented in: 
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

%% Process inputs
% Tr1 = Tr1;% - mean(Tr1(1:20));
% Tr2 = Tr2;% - mean(Tr2(1:20));
M = range(Tr2)/range(Tr1);
Tr1 = Tr1*M;


tof = tof(:);
Tr1 = Tr1(:);
Tr2 = Tr2(:);
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
optimiseHP = 'none';  % if true optimises the lengthscale for each Transmission spectra
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

%% fit exponentials to right and left of both Tr1 and Tr2

%% Proj
%% 1) fit to the far right of the edge where B = 1, so only fit exp([-(a0+b0.*t)])
fit1Right = @(p,x) exp(-(p(1) + p(2).*x));
[p,~,~,~,~,~,~] = lsqcurvefit(fit1Right,[a00;b00],tof(opts.endIdx(1):opts.endIdx(2)),Tr1(opts.endIdx(1):opts.endIdx(2)),[],[],optionsFit);
a0 = p(1); b0 = p(2);
%% 2) fit to the far left of the edge where B = 0;
fit1Left = @(p,x) exp(-(a0 + b0.*x)).*exp(-(p(1)+p(2).*x));
[p,~,~,~,~,~,~] = lsqcurvefit(fit1Left,[a_hkl0;b_hkl0],tof(opts.startIdx(1):opts.startIdx(2)),Tr1(opts.startIdx(1):opts.startIdx(2)),[],[],optionsFit);
a_hkl = p(1); b_hkl = p(2);

g1Proj = @(x) exp(-(a0 + b0.*x)).*exp(-(a_hkl+b_hkl.*x));
g2Proj = @(x) exp(-(a0 + b0.*x));

%----------------------------------------------------------------------------------%
if ~isfield(opts,'gD0')
%% D0
%% 1) fit to the far right of the edge where B = 1, so only fit exp([-(a0+b0.*t)])
fit1Right = @(p,x) exp(-(p(1) + p(2).*x));
[p,~,~,~,~,~,~] = lsqcurvefit(fit1Right,[a00;b00],tof(opts.endIdx(1):opts.endIdx(2)),Tr2(opts.endIdx(1):opts.endIdx(2)),[],[],optionsFit);
a0 = p(1); b0 = p(2);
%% 2) fit to the far left of the edge where B = 0;
fit1Left = @(p,x) exp(-(a0 + b0.*x)).*exp(-(p(1)+p(2).*x));
[p,~,~,~,~,~,~] = lsqcurvefit(fit1Left,[a_hkl0;b_hkl0],tof(opts.startIdx(1):opts.startIdx(2)),Tr2(opts.startIdx(1):opts.startIdx(2)),[],[],optionsFit);
a_hkl = p(1); b_hkl = p(2);

g1D0 = @(x) exp(-(a0 + b0.*x)).*exp(-(a_hkl+b_hkl.*x));
g2D0 = @(x) exp(-(a0 + b0.*x));

%% D0
yD0 = (Tr2-g1D0(tof));
x=tof;
sig_m_D0 = std([Tr2(opts.endIdx(1):opts.endIdx(2)) - g2D0(tof(opts.endIdx(1):opts.endIdx(2)));...
    Tr2(opts.startIdx(1):opts.startIdx(2))-g1D0(tof(opts.startIdx(1):opts.startIdx(2)))]);
end
xt = linspace(tof(opts.startIdx(2)),tof(opts.endIdx(1)),nx)';

%----------------------------------------------------------------------------------%
%% fit Transition as GP
%% calculate sig_m for both Tr1 and Tr2
%% PROJ
yProj = (Tr1-g1Proj(tof));
x=tof;
sig_m_Proj = std([Tr1(opts.endIdx(1):opts.endIdx(2)) - g2Proj(tof(opts.endIdx(1):opts.endIdx(2)));...
    Tr1(opts.startIdx(1):opts.startIdx(2))-g1Proj(tof(opts.startIdx(1):opts.startIdx(2)))]);

%%
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
PhiProj = Phi .* (g2Proj(x) - g1Proj(x));
% PhiD0 = Phi .* (g2D0(x) - g1D0(x));
[~,m]= size(PhiProj);

GammaProj = [PhiProj./sig_m_Proj;diag(1./sqrt(SLambda))];
% GammaD0 = [PhiD0./sig_m_D0;diag(1./sqrt(SLambda))];

RProj = triu(qr(GammaProj));
% RD0   = triu(qr(GammaD0));

CZProj = RProj(1:m,1:m);

if ~isfield(opts,'gD0')
PhiD0 = Phi .* (g2D0(x) - g1D0(x));
GammaD0 = [PhiD0./sig_m_D0;diag(1./sqrt(SLambda))];
RD0   = triu(qr(GammaD0));
CZD0 = RD0(1:m,1:m);
LSoptsT.TRANSA = true; LSoptsT.UT = true; LSopts.TRANSA = false; LSopts.UT = true;
v_basisD0 = (linsolve(CZD0,linsolve(CZD0,PhiD0'*(yD0./sig_m_D0.^2),LSoptsT),LSopts));
festpD0     = Phi_p*v_basisD0;  %
opts.gD0 = dPhi_T * v_basisD0;
end
gD0 = opts.gD0;

LSoptsT.TRANSA = true; LSoptsT.UT = true; LSopts.TRANSA = false; LSopts.UT = true;


v_basisProj = (linsolve(CZProj,linsolve(CZProj,PhiProj'*(yProj./sig_m_Proj.^2),LSoptsT),LSopts));
% v_basisD0 = (linsolve(CZD0,linsolve(CZD0,PhiD0'*(yD0./sig_m_D0.^2),LSoptsT),LSopts));


festpProj   = Phi_p*v_basisProj;  %
% festpD0     = Phi_p*v_basisD0;  %


gProj = dPhi_T * v_basisProj;
% gD0 = dPhi_T * v_basisD0;
% alternatively
gamma = linsolve(CZProj,eye(m),LSoptsT);
SigV = triu(qr(gamma));
svProj = v_basisProj  + SigV.'*randn(m,ns);      % sample the coefficients
sg  = dPhi_T * svProj;
%% calculate normalised cross correlation for all samples
XC = nan(nx*2-1,ns);
for i =1:ns
    [XC(:,i),lags] = xcorr(sg(:,i),gD0,[],'normalized');
%     [XC(:,i),lags] = xcorr(sg(:,i),sg2(:,i),[],'normalized');
end

dt = xt(2)-xt(1);
%% take max
[~,I] = max(XC,[],1);
%% calculate sample mean and standard deviation
d_all = lags(I) * dt;

deltaD = mean(d_all);
sigma  = std(d_all);
TrFit  = festpProj;

fitinfo.lengthscale = l;
fitinfo.d_all = d_all;
end