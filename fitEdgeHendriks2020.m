function [edgePos,sigma,TrFit] = fitEdgeHendriks2020(Tr,tof,opts)
%fitEdgeHendriks2020 fits a bragg-edge using the method presented in:
%   TODO: Insert arxiv link when paper is written.
%
% Inputs:
%   - Tr is a Nx1 double containing the normalised transmisssion curve
%   for a single projection
%   - tof is an Nx1 array of wave-lengths or time-of-flight.
%   - options is a structure containing
%       opts.a00    :   Initial guess
%       opts.b00    :   Initial guess
%       opts.a_hkl  :   Initial guess
%       opts.b_hkl  :   Initial guess
%
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
% Last modified: 08/01/2020
% This program is licensed under GNU GPLv3, see LICENSE for more details.

%% least squares fitting options
optionsFit              = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');
optionsFit.Algorithm    = 'Levenberg-Marquardt';
optionsFit.Jacobian     = 'off';
optionsFit.Display      = 'off';
%% Initial guess
a00 = 0.5;
b00 = 0.5;
a_hkl0 = 0.5;
b_hkl0 = 0.5;
% p00 = [0.0187,0.006,0.008,1];
p00 = [mean([opts.startRange(2) opts.endRange(1)]),... % Edge location
    (tof(2)-tof(1))*1e3,... % width
    (tof(2)-tof(1))*1e3,... % assymetry
    0,...   %pedistool
    0.5];     %slope

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
%TODO add hyperparameters to struct
%TODO add number of samples to struct
%TODO move to opts struct
ns = 3000;
nx = 2500;
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
sig_f = 1;
% l = 0.02;
l = 9e-5;

%GP
xt = linspace(tof(opts.startIdx(2)),tof(opts.endIdx(1)),nx)';

x = tof.';
K = sig_f^2 * exp(-0.5*(x - x').^2/l^2) .* ((g2(x) - g1(x)) .*(g2(x') - g1(x')));
Kyy = K + eye(ny)*sig_m^2;

Kfy = sig_f^2 * exp(-0.5*(xt - x').^2/l^2).*(g2(x') - g1(x'));
dKfy = -(xt - x')/l^2 .* Kfy;
ddKff = sig_f^2 * (1 - (xt-xt').^2/l^2)/l^2 .* exp(-0.5*(xt - xt').^2/l^2);

% estimated at the original wavelengths for plotting purposes
Kfyp = sig_f^2 * exp(-0.5*(x - x').^2/l^2).*(g2(x') - g1(x'));
festp = Kfyp*(Kyy\y);

C = chol(Kyy,'upper');
g = dKfy*(C\(C'\y));

alpha = (C'\dKfy');
V = ddKff- alpha'*alpha;



[sV, p] = chol(V+1e-10*eye(size(V)),'upper');
if p
    warning('was not about to use chol')
    sV = sqrtm(V+1e-10*eye(size(V)));
    sg = g + sV*randn(nx,ns);     % notice this one is not tranposed (for good reason)
else
    sg = g + sV.'*randn(nx,ns);
end

[~,Is] = max([g sg]);

sLams = xt(Is);
%% Collect Results
edgePos = mean(sLams);
sigma = std(sLams);
TrFit = festp;
% Calculate d-spacings and confidence
% d_cell{k}(i)=mean(sLams);
% std_cell{k}(i)=std(sLams);  %% not sure how to do this?
% edge_shape{k}(i,:) = festp;

end