function [edgePos,sigma,TrFit,fitinfo] = fitEdge5ParamMethod(Tr,tof,opts)
%fitEdge5ParamMethod fits a bragg-edge using the method used by: Tremsin,
%   A. S., Gao, Y., Dial, L. C., Grazzi, F., Shinohara, T., 2016.
%   Investigation of microstructure in additive manufactured inconel 625 by
%   spatially resolved neutron transmission spectroscopy. Science and
%   Technology of advanced MaTerialS 17 (1), 324?336.
%   https://doi.org/10.1080/14686996.2016.1190261
%
% Inputs:
%   - Tr is a Nx1 double containing the normalised transmisssion curve
%   for a single projection
%   - tof is an Nx1 array of wave-lengths or time-of-flight.
%   - options is a structure containing
%       opts.t_hkl0 :   Initial guess for edge location
%       opts.sigma0 :   Initial guess for gaussian broadening term
%       opts.tau0   :   Initial guess for exponential decay term
%       opts.C10    :   Initial guess for pedistool
%       opts.C20    :   Initial guess for slope
%
% Outputs:
%   - edgePos is the location of the braggEdge
%   - sigma is the estimated standard deviation
%   - TrFit is is the Bragg edge model evaluated at tof
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
% Last modified: 15/01/2020
% This program is licensed under GNU GPLv3, see LICENSE for more details.

%% least squares fitting options
optionsFit              = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');
optionsFit.Algorithm    = 'Levenberg-Marquardt';
optionsFit.Jacobian     = 'off';
optionsFit.Display      = 'off';
%% Initial guess

p00 = [mean(opts.range),... % Edge location
    1e-5,... % width
    1e-5,... % assymetry 
    min(Tr),...   %pedistool
    0.5*(max(Tr)-min(Tr))];     %slope

if isfield(opts,'t_hkl0')
    p00(1) = opts.t_hkl0;
end
if isfield(opts,'sigma0')
    p00(2) = opts.sigma0;
end
if isfield(opts,'tau0')
    p00(3) = opts.tau0;
end
if isfield(opts,'C10')
    p00(4) = opts.C10;
end
if isfield(opts,'C20')
    p00(5) = opts.C20;
end
%% Fit edge
idx = opts.rangeIdx(1):opts.rangeIdx(2);
fitMe = @(p,x) edgeModel(p,x);
[p,resnorm,residual,~,~,~,J] = lsqcurvefit(fitMe,p00,tof(idx),Tr(idx),[],[],optionsFit);
%% Collect Results
edgePos = p(1);
ci = nlparci(p,residual,'jacobian',J); % confidence intervals
sigma = (ci(1,2)-ci(1,1))/4;
TrFit = fitMe(p,tof);

fitinfo.resnorm = resnorm;
fitinfo.edgewidth = p(2);
fitinfo.egdgeassymetry = p(3);
end

function [edge_spect] = edgeModel(params,t)
t_hkl = params(1);      % edge location
sigma = params(2);      % width (broadening)
tau = params(3);        % assymetry
C1  = params(4);        % height
C2  = params(5);        % slope

%% As presented in : Tremsin A. S., Gao, Y., Dial, L. C., Grazzi, F., Shinohara, T., 2016.
%   Investigation of microstructure in additive manufactured inconel 625 by
%   spatially resolved neutron transmission spectroscopy. Science and
%   Technology of advanced MaTerialS 17 (1), 324?336.
%   https://doi.org/10.1080/14686996.2016.1190261
%
% edge_spect = C1+C2.*(erfc(-(t-t_hkl)./(sqrt(2)*sigma))...
%     - exp(-(t-t_hkl)./tau + sigma^2./(2*tau.^2))...
%     .*erfc(-(t-t_hkl)./(sqrt(2)*sigma)+sigma./(sqrt(2)*tau)));
% 
% The 1/sqrt(2) in the final term was confirmed as a mis-print by Tremsin
% the intended equation is:
edge_spect = C1+C2.*(erfc(-(t-t_hkl)./(sqrt(2)*sigma))...
    - exp(-(t-t_hkl)./tau + sigma^2./(2*tau.^2))...
    .*erfc(-(t-t_hkl)./(sqrt(2)*sigma)+sigma./tau)); 
end