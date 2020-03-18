function [edgePos,sigma,TrFit] = fitEdgeGPMethod(Tr,tof,opts)
%fitEdgeGPMethod fits a bragg-edge using the method presented in:
%   TODO: Insert Bib entry
%   TODO: Insert arxiv link when paper is written.
%
% Inputs:
%   - Tr is a Nx1 double containing the normalised transmisssion curve
%   for a single projection
%   - tof is an Nx1 array of wave-lengths or time-of-flight.
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
%
%
% Copyright (C) 2020 The University of Newcastle, Australia
% Authors:
%   Nicholas O'Dell <Nicholas.Odell@newcastle.edu.au>
%   Johannes Hendriks <Johannes.Hendriks@newcastle.edu.au>
% Last modified: 18/03/2020
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
useInterp = true;   % uses an interpolation procedure to reduce the size of matrix inversion

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
if isfield(opts,'l')
    l = opts.l;
end
if isfield(opts,'ns')
    ns = opts.ns;
end
if isfield(opts,'nx')
    nx = opts.nx;
end
if isfield(opts,'useInterp')
   useInterp = opts.useInterp; 
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
ny = length(tof);
sig_m = std([Tr(opts.endIdx(1):opts.endIdx(2)) - g2(tof(opts.endIdx(1):opts.endIdx(2))),...
    Tr(opts.startIdx(1):opts.startIdx(2))-g1(tof(opts.startIdx(1):opts.startIdx(2)))]);
% Hyperparameters

%GP
if useInterp
   nh = nx;
   nx = length(tof); 
   xt_interp = linspace(tof(opts.startIdx(2)),tof(opts.endIdx(1)),nh)';
end

xt = linspace(tof(opts.startIdx(2)),tof(opts.endIdx(1)),nx)';

x = tof.';
K = sig_f^2 * exp(-0.5*(x - x').^2/l^2) .* ((g2(x) - g1(x)) .*(g2(x') - g1(x')));
Kyy = K + eye(ny)*sig_m^2;

Kfy = sig_f^2 * exp(-0.5*(xt - x').^2/l^2).*(g2(x') - g1(x'));
dKfy = -(xt - x')/l^2 .* Kfy;
ddKff = sig_f^2 * (1 - (xt-xt').^2/l^2)/l^2 .* exp(-0.5*(xt - xt').^2/l^2);
Kfyp = sig_f^2 * exp(-0.5*(x - x').^2/l^2).*(g2(x') - g1(x'));
 
C = chol(Kyy,'upper');

festp = Kfyp*(C\(C.'\y));           % estimated edge shape

g = dKfy*(C\(C.'\y));

alpha = (C.'\dKfy');
V = ddKff- alpha'*alpha;



[sV, p] = chol(V+1e-10*eye(size(V)),'upper');
if p
    warning('was not about to use chol')
    sV = sqrtm(V+1e-10*eye(size(V)));
    sg = g + sV*randn(nx,ns);     % notice this one is not tranposed (for good reason)
else
    sg = g + sV.'*randn(nx,ns);
end

if useInterp
    g_interp = interp1(xt,g,xt_interp,'v5cubic');
    sg_interp = interp1(xt,sg,xt_interp,'v5cubic');
    [~,Is] = max([g_interp sg_interp]);
    sLams = xt_interp(Is);
else
    [~,Is] = max([g sg]);
    sLams = xt(Is);
end


%% Collect Results
edgePos = mean(sLams);
sigma = std(sLams);
TrFit = exp(-a0-b0*tof).*...
	(exp(-a_hkl-b_hkl*tof) + (1-exp(-a_hkl - b_hkl*tof)) .*festp.');

end