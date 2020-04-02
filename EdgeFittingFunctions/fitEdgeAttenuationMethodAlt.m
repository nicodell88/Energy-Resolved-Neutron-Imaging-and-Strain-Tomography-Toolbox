function [edgePos,sigma,TrFit,fitinfo] = fitEdgeAttenuationMethodAlt(Tr,tof,opts)
%fitEdgeAttenuationMethod fits a bragg-edge using the method presented in:
%   Santisteban, J., Edwards, L., Steuwer, A., Withers, P., 2001.
%   Time-of-flight neutron transmission diffraction. Journal of applied
%   crystallography 34 (3), 289?297.
%   https://onlinelibrary.wiley.com/doi/pdf/10.1107/S0021889801003260
%   However, instead of using the edge model presented in Santisteban's
%   work, an alternative edge model is used. This edge model is given on
%   page 34 of
%   Vogel, Sven. A Rietveld-approach for the analysis of neutron time-of-flight
%   transmission data. Diss. Christian-Albrechts Universität Kiel, 2000.
%   This edge model is given by the integral of a convolution between a
%   Gaussian and two back to back exponentials with different rise and
%   decay times alpha and beta
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
%       opts.t_hkl0 :   Initial guess for edge location
%       opts.sigma0 :   Initial guess for gaussian broadening term
%       opts.alpha  :   Initial guess for exponential rise term
%       opts.beta   :   Initial guess for exponential decay term
%
% Outputs:
%   - edgePos is the location of the braggEdge
%   - sigma is the estimated standard deviation
%   - TrFit is the Bragg edge model evaluated at tof
%   - fitinfo contains additional information about the fit
%       fitinfo.resnorm     : the squared 2 norm of the residual as
%                               calcualted by lsqcurvefit
%       fitinfo.edgewidth   : edge width parameter from the attenuation model 
%       fitinfo.alpha       : exponential rise time
%       fitinfo.beta        : exponential decay time
%
%
% Copyright (C) 2020 The University of Newcastle, Australia
% Authors:
%   Johannes Hendriks <Johannes.hendriks@newcastle.edu.au>
% Last modified: 25/03/2020
% This program is licensed under GNU GPLv3, see LICENSE for more details.

%% least squares fitting options
optionsFit              = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');
optionsFit.Algorithm    = 'Levenberg-Marquardt';
optionsFit.Jacobian     = 'off';
optionsFit.Display      = 'off';
optionsFit.FunctionTolerance = 1e-10;
optionsFit.StepTolerance = 1e-10;
%% Initial guess
a00 = 0.5;
b00 = 0.5;
a_hkl0 = 0.5;
b_hkl0 = 0.5;
% p00 = [0.0187,0.006,0.008,1];
p00 = [mean([opts.startRange(2) opts.endRange(1)]),... % Edge location
    (tof(2)-tof(1))*1e1,... % width
    1/((tof(2)-tof(1))*4),...% alpha0
    1/((tof(2)-tof(1)))];  %beta0

v00 = [mean([opts.startRange(2) opts.endRange(1)]),... % Edge location
    (tof(2)-tof(1))*1e3,... % width
    (tof(2)-tof(1))*1e3]; % assymetry 

if isfield(opts,'alpha0')
    a00 = opts.a00;
end
if isfield(opts,'beta0')
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
    v00(1) = opts.t_hkl0;
end
if isfield(opts,'sigma0')
    p00(2) = opts.sigma0;
    v00(2) = opts.sigma0;
end
if isfield(opts,'alpha0')
    p00(3) = opts.alpha0;
end
if isfield(opts,'beta0')
    p00(4) = opts.beta0;
end
%% Fit edge
try
% 1) fit to the far right of the edge where B = 1, so only fit exp([-(a0+b0.*t)])
fit1 = @(p,x) exp(-(p(1) + p(2).*x));
[p,~,~,~,~,~,~] = lsqcurvefit(fit1,[a00;b00],tof(opts.endIdx(1):opts.endIdx(2)),Tr(opts.endIdx(1):opts.endIdx(2)),[],[],optionsFit);
a0 = p(1); b0 = p(2);
% 2) fit to the far left of the edge where B = 0;
fit2 = @(p,x) exp(-(a0 + b0.*x)).*exp(-(p(1)+p(2).*x));
[p,~,~,~,~,~,~] = lsqcurvefit(fit2,[a_hkl0;b_hkl0],tof(opts.startIdx(1):opts.startIdx(2)),Tr(opts.startIdx(1):opts.startIdx(2)),[],[],optionsFit);
a_hkl = p(1); b_hkl = p(2);
% 3) now fit the area around the edge

% find some initial values for t_hkl and sigma
fit3 = @(p,x) edgeModel([p a0 b0 a_hkl b_hkl],x);
[p] = lsqcurvefit(fit3,v00,tof,Tr,[],[],optionsFit);
sig = p(2);

% get initial values for alpha and beta while holding sigma constant
fit4 = @(p,x) edgeModelAlt([p(1) sig p(2:3) a0 b0 a_hkl b_hkl],x);
[p,~,~,~,~,~,~] = lsqcurvefit(fit4,[p00(1) p00(3:4)],tof(opts.startIdx(2):opts.endIdx(1)),Tr(opts.startIdx(2):opts.endIdx(1)),[],[],optionsFit);

% now let all be free
fit4 = @(p,x) edgeModelAlt([p a0 b0 a_hkl b_hkl],x);
p00 = [p(1) sig p(2:3)];
[p,~,residual,~,~,~,J] = lsqcurvefit(fit4,p00,tof(opts.startIdx(2):opts.endIdx(1)),Tr(opts.startIdx(2):opts.endIdx(1)),[],[],optionsFit);
catch
    edgePos = NaN;
    sigma = NaN;
    TrFit = nan(size(tof));
    fitinfo.resnorm = NaN;
    fitinfo.edgewidth = NaN;
    fitinfo.alpha = NaN;
    fitinfo.beta = NaN;
    return  
end
%% Collect Results
edgePos = p(1);
ci = nlparci(p,residual,'jacobian',J); % confidence intervals
sigma = (ci(1,2)-ci(1,1))/4;
TrFit = fit4(p,tof);
TrFit(TrFit<-0.99) = NaN;
TrFit(TrFit>1.99) = NaN;

inds = ~isnan(TrFit);
residual = Tr(inds) - TrFit(inds);

% figure(1)
% clf
% plot(Tr)
% hold on
% plot(TrFit)
% hold off

fitinfo.resnorm = sum(residual.^2);
fitinfo.edgewidth = p(2);
fitinfo.alpha = p(3);
fitinfo.beta = p(4);
end

function [edge_spect] = edgeModel(params,t)
t_hkl = params(1);      % edge location
sigma = params(2);      % width (broadening)
tau = params(3);        % assymetry
% v = params(4);

a0 = params(4);
b0 = params(5);
a_hkl = params(6);
b_hkl = params(7);

B = 1/2.*(erfc(-(t-t_hkl)./(sqrt(2)*sigma))...
    - exp(-(t-t_hkl)./tau + sigma^2./(2*tau.^2))...
    .*erfc(-(t-t_hkl)./(sqrt(2)*sigma)+sigma./tau));

A1 = exp(-(a0+b0.*t));
A2 = exp(-(a_hkl+b_hkl.*t));

edge_spect = A1.*( A2 +(1-A2).*(B));

end

function [edge_spect] = edgeModelAlt(params,t)
t_hkl = params(1);      % edge location
sigma = params(2);      % width (broadening)
alpha = params(3);        % exponential rise rate
beta = params(4);        % exponential decay rate

a0 = params(5);
b0 = params(6);
a_hkl = params(7);
b_hkl = params(8);

delta = (t_hkl-t);
w = delta/sqrt(2)/sigma;
u = alpha/2*(alpha*sigma^2+2*delta);
v = beta/2*(beta*sigma^2-2*delta);
y = (alpha*sigma^2+delta)/sqrt(2)/sigma;
z = (beta*sigma^2-delta)/sqrt(2)/sigma;


B = 0.5*erfc(w) - (beta*exp(u).*erfc(y) - alpha*exp(v).*erfc(z))./(2*(alpha+beta));



A1 = exp(-(a0+b0.*t));
A2 = exp(-(a_hkl+b_hkl.*t));

edge_spect = A1.*( A2 +(1-A2).*(B));

edge_spect(isnan(edge_spect)) = -1.0;
edge_spect(isinf(edge_spect)) = 2.0;
if any(isnan(edge_spect)) || any(isinf(edge_spect))
    disp('bad')
end

end