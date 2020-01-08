function [edgePos,sigma,TrFit] = fitEdgeTremsin2011(Tr,tof,opts)
%fitEdgeTremsin2011 fits a bragg-edge using the method presented in:
%   https://doi.org/10.1080/14686996.2016.1190261
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
%       opts.p00    :   Initial guess
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
if isfield(opts,'t_hkl0')
    p00(1) = opts.t_hkl0;
end
if isfield(opts,'sigma')
    p00(2) = opts.sigma;
end
if isfield(opts,'tau')
    p00(3) = opts.tau;
end
if isfield(opts,'C1')
    p00(4) = opts.C1;
end
if isfield(opts,'C2')
    p00(5) = opts.C2;
end
%% Fit edge
% 1) fit to the far right of the edge where B = 1, so only fit exp([-(a0+b0.*t)])
fit1 = @(p,x) exp(-(p(1) + p(2).*x));
[p,~,~,~,~,~,~] = lsqcurvefit(fit1,[a00;b00],tof(opts.endIdx(1):opts.endIdx(2)),Tr(opts.endIdx(1):opts.endIdx(2)),[],[],optionsFit);
a0 = p(1); b0 = p(2);
% 2) fit to the far left of the edge where B = 0;
fit2 = @(p,x) exp(-(a0 + b0.*x)).*exp(-(p(1)+p(2).*x));
[p,~,~,~,~,~,~] = lsqcurvefit(fit2,[a_hkl0;b_hkl0],tof(opts.startIdx(1):opts.startIdx(2)),Tr(opts.startIdx(1):opts.startIdx(2)),[],[],optionsFit);
a_hkl = p(1); b_hkl = p(2);
% 3) now fit the area around the edge
fit3 = @(p,x) edgeModel([p a0 b0 a_hkl b_hkl],x);
[p,~,residual,~,~,~,J] = lsqcurvefit(fit3,p00,tof,Tr,[],[],optionsFit);
%% Collect Results
edgePos = p(1);
ci = nlparci(p,residual,'jacobian',J); % confidence intervals
sigma = (ci(1,2)-ci(1,1))/4;
TrFit = fit3(p,tof);
end

function [edge_spect] = edgeModel(params,t)
t_hkl = params(1);      % edge location
sigma = params(2);      % width (broadening)
tau = params(3);        % assymetry
C1  = params(4);        % height
C2  = params(5);        % slope

a0 = params(6);
b0 = params(7);
a_hkl = params(8);
b_hkl = params(9);

B = C1+C2.*(erfc(-(t-t_hkl)./(sqrt(2)*sigma))...
    - exp(-(t-t_hkl)./tau + sigma^2./(2*tau.^2))...
    .*erfc(-(t-t_hkl)./(sqrt(2)*sigma)+sigma./tau));

A1 = exp(-(a0+b0.*t));
A2 = exp(-(a_hkl+b_hkl.*t));

edge_spect = A1.*( A2 +(1-A2).*(B));
end