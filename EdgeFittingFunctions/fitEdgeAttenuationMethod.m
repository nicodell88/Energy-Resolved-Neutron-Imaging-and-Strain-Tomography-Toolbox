function [edgePos,sigma,TrFit,fitinfo] = fitEdgeAttenuationMethod(Tr,tof,opts)
%fitEdgeAttenuationMethod fits a bragg-edge using the method presented in:
%   Santisteban, J., Edwards, L., Steuwer, A., Withers, P., 2001.
%   Time-of-flight neutron transmission diffraction. Journal of applied
%   crystallography 34 (3), 289?297.
%   https://onlinelibrary.wiley.com/doi/pdf/10.1107/S0021889801003260
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
%       opts.tau0   :   Initial guess for exponential decay term
%
% Outputs:
%   - edgePos is the location of the braggEdge
%   - sigma is the estimated standard deviation
%   - TrFit is the Bragg edge model evaluated at tof
%   - fitinfo contains additional information about the fit
%       fitinfo.resnorm     : the squared 2 norm of the residual as
%                               calcualted by lsqcurvefit
%       fitinfo.edgewidth   : edge width parameter from the attenuation model 
%       fitinfo.egdgeassymetry : edge assymetry parameter of atten model
%
%
% Copyright (C) 2020 The University of Newcastle, Australia
% Authors:
%   Nicholas O'Dell <Nicholas.Odell@newcastle.edu.au>
%   Johannes Hendriks <Johannes.hendriks@newcastle.edu.au>
% Last modified: 10/01/2020
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
    (tof(2)-tof(1))*1e3]; % assymetry 

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
if isfield(opts,'sigma0')
    p00(2) = opts.sigma0;
end
if isfield(opts,'tau0')
    p00(3) = opts.tau0;
end

try
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
[p,resnorm,residual,~,~,~,J] = lsqcurvefit(fit3,p00,tof,Tr,[],[],optionsFit);
catch
    edgePos = NaN;
    sigma = NaN;
    TrFit = nan(size(tof));
    fitinfo.resnorm = NaN;
    fitinfo.edgewidth = NaN;
    fitinfo.alpha = NaN;
    fitinfo.beta = NaN;
end
%% Collect Results
edgePos = p(1);
ci = nlparci(p,residual,'jacobian',J); % confidence intervals
sigma = (ci(1,2)-ci(1,1))/4;
TrFit = fit3(p,tof);

fitinfo.resnorm = resnorm;
fitinfo.edgewidth = p(2);
fitinfo.egdgeassymetry = p(3);
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