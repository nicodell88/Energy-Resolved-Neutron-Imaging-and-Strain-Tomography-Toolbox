function [deltaD,sigma,TrFit,fitinfo,opts] = crossCorrMethod(Tr1,Tr2,tof,opts)
%crossCorrMethod Calculates the "shift" in wavelength or time of flight
%	between two sets of transmission data. Through an unstressed sample and a
%	stressed sample. This method is presented in: Ramadhan, R.S., Kockelmann,
%	W., Minniti, T., Chen, B., Parfitt, D., Fitzpatrick, M.E. and Tremsin,
%	A.S., 2019. Characterization and application of Bragg-edge transmission
%	imaging for strain measurement and crystallographic analysis on the IMAT
%	beamline. Journal of Applied Crystallography, 52(2).
%	https://journals.iucr.org/j/issues/2019/02/00/ks5612/ks5612.pdf
%
% Inputs:
%   - Tr1 is a Nx1 double containing raw transmisssion curve
%   for a single projection (stresseed)
% 	- Tr2 is a Nx1 double containing the raw transmisssion curve
%   for a single projection (unstressed)
%   - tof is an Nx1 array of wave-lengths or time-of-flight.
%   - options is a structure containing
%       opts.range      :   a 2 element vector contianing the window for
%                           fitting [start end] in the same units as tof.
%       opts.order      :   order for SGOLAY filter
%       opts.frame      :   framelength for SGOLAY filter
%       opts.peakWindow :   sample window for peak fitting voigt function.
%
% Outputs:
%	- deltaD is shift between the Bragg-Edge in the un-stressed and
%       stressed sample.
%   - sigma is the estimated standard deviation
%
% Copyright (C) 2020 The University of Newcastle, Australia
% Authors:
%   Nicholas O'Dell <Nicholas.Odell@newcastle.edu.au>
%   Johannes Hendriks <Johannes.Hendriks@newcastle.edu.au>
% Last modified: 20/03/2020
% This program is licensed under GNU GPLv3, see LICENSE for more details.

Tr1 = Tr1 - mean(Tr1(1:20));
Tr2 = Tr2 - mean(Tr2(1:20));
M = mean(Tr2(end-20:end))/mean(Tr1(end-20:end));
Tr1 = Tr1*M;

Tr1 = Tr1(:);
Tr2 = Tr2(:);
tof = tof(:);
%% Curve fitting options for fitting voigt function
optionsFit              = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');
optionsFit.Algorithm    = 'Levenberg-Marquardt';
optionsFit.Jacobian     = 'off';
optionsFit.Display      = 'off';
optionsFit.OptimalityTolerance = 1e-9;

%% populate opts structure
if nargin>3
    if ~isfield(opts,'range')
        opts.range = [tof(1) tof(end)];
    end
    
    if ~isfield(opts,'order')
        opts.order = 3;
    end
    
    if ~isfield(opts,'frame')
        opts.frame = 35;
    end
    if ~isfield(opts,'peakWindow')
        opts.peakWindow = 20;
    end
    
    if ~isfield(opts,'y0')
        opts.y0 = 0;
    end
    if ~isfield(opts,'xc')
        opts.xc = 0;
    end
    if ~isfield(opts,'A')
        opts.A = exp(4.5);
    end
    if ~isfield(opts,'Wg')
        opts.Wg = exp(4.3);
    end
    if ~isfield(opts,'Wl')
        opts.Wl = exp(4.38);
    end
    if ~isfield(opts,'Mu')
        opts.Mu = exp(-3.5);
    end
    opts.p00 = [...
        opts.y0
        opts.xc
        log(opts.A)
        log(opts.Wg)
        log(opts.Wl)
        log(opts.Mu)];
    
else
    opts.range = [tof(1) tof(end)];
    opts.order = 5;
    opts.frame = 7;
    opts.peakWindow = 20;
    
    opts.p00 = [  ...
        15.2407     %y0
        0.0         %xc
        0      %log(A)
        -4.1025     %log(Wg)
        -3.5159     %log(Wl)
        0.0876];    %log(Mu)
end


%% Inspect data
[~,opts.fitIdx] = min((tof(:).'-opts.range(:)).^2,[],2);
if any(abs(diff(diff(tof(opts.fitIdx(1):opts.fitIdx(2))))) > 1e-5)
    warning('Data must be uniformly sampled.')
end

%% Smooth data
[~,g] = sgolay(opts.order,opts.frame);

p = 1;  %1st derivative
dt = tof(opts.fitIdx(2))-tof(opts.fitIdx(2)-1);
dTr1 = conv(Tr1(opts.fitIdx(1):opts.fitIdx(2)),factorial(p)/(-dt)^p * g(:,p+1));
dTr2 = conv(Tr2(opts.fitIdx(1):opts.fitIdx(2)),factorial(p)/(-dt)^p * g(:,p+1));

dTr1 = dTr1(opts.frame:(end-opts.frame));
dTr2 = dTr2(opts.frame:(end-opts.frame));

%% Calculate cross correlation
[C,Lags] = xcorr(dTr1,dTr2,[],'normalized');
% X = Lags*dt*500;
X = Lags*dt;

TrFit = nan(size(tof));

%% Fit peak
window = floor(opts.peakWindow/2);
[~,idx] = max(C);
idxFit = (-window:1:window) + idx;
idxFit(idxFit > length(C)) = [];
idxFit(idxFit < 1) = [];

opts.p00(2) = X(idx);

fit1 = @(p,x) pseudoVoigt([p],x);
[p,~,residual,~,~,~,J] = lsqcurvefit(fit1,opts.p00,X(idxFit),C(idxFit),[],[],optionsFit);
ci = nlparci(p,residual,'jacobian',J); % confidence intervals
%%
% p = [...
%     0
%     -50
%     2.6
%     3
%     2
%     0]

figure(2)
clf
plot(X(idxFit),C(idxFit))
hold on
x = X(idxFit);
y = fit1(p,x);
plot(x,y)

% pause
%% Extract results
% deltaD = p(2)/500;
deltaD = p(2);
sigma = (ci(2,2)-ci(2,1))/4/500;
fitinfo.resnorm = residual;
fitinfo.p = p;
% fitinfo.edgewidth = p(2);
% fitinfo.egdgeassymetry = p(3);
end

function y = pseudoVoigt(p,x)

y0  = p(1);
xc  = p(2);
A   = exp(p(3));
Wg  = exp(p(4));
Wl  = exp(p(5));
Mu  = exp(p(6));

y = y0 + A*(Mu * 2/pi * Wl./(4*(x-xc).^2 + Wl^2) +...
    (1-Mu)*(sqrt(4*log(2))./(sqrt(pi)*Wg))*exp(-(4*log(2)*(x-xc).^2)./(Wg^2)) );

y=y(:);
end

