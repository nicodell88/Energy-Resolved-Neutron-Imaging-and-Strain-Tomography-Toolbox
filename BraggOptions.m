function opts = BraggOptions(tof,method,varargin)
%BraggOptions(tof,method,...) generates an options structure which contains the
%neccecary options for edge fitting. This structure is used during any task
%requiring edge fitting.
%   opts = BraggOptions(tof) specifies 'gp' as the default edge fitting
%   method and generates the default options structure. tof is a vector of
%   time of flight spectra.
%
%   opts = BraggOptions(tof,'attenuation') generates the default
%   structure for the Sanitisteban edge fitting method. Further options can
%   be specified using name-value pairs. E.g., the default edge location
%   can be specified using: 
%   opts = BraggOptions(tof,'attenuation','t_hkl0', 4.17)
%
%   More information on the name-value pairs for each method can be found
%   in the help documentation for each edge fitting function.
%
%   See also fitEdges, fitEdgeAttenuationMethod,
%   fitEdgeAttenuationMethodAlt, fitEdge5ParamMethod, fitEdgeGPMethod,
%   fitEdgeGPCCMethod, crossCorrMethod.
%

% Copyright (C) 2020 The University of Newcastle, Australia
% Authors:
%   Nicholas O'Dell <Nicholas.Odell@newcastle.edu.au>
% Last modified: 23/06/2020
% This program is licensed under GNU GPLv3, see LICENSE for more details.

validateattributes((tof),{'numeric'},{'increasing'})

if~exist('method','var') || isempty(method)
    method = 'gp';
end

p = inputParser;
%% Method
defaultMethod = 'gp';
expectedMethods = {'attenuation','attenuationalt','5param','crosscorr','gp','gpcc','gpcc2'};
addOptional(p,'method',defaultMethod, ...
    @(x) any(validatestring(x,expectedMethods)));
%% Plot
default = false;
addParameter(p,'plot',default, ...
    @(x) islogical(x));
%% Parallel
default = false;
addParameter(p,'Par',default, ...
    @(x) islogical(x));
%% Parse
parse(p,method);
%% Add options specific to each method
%% a00
default = 0.5;
addParameter(p,'a00',default, ...
    @(x) validateattributes((x),{'numeric'},{'scalar'}));
%% b00
default = 0.5;
addParameter(p,'b00',default, ...
    @(x) validateattributes((x),{'numeric'},{'scalar'}));
%% a_hkl0
default = 0.5;
addParameter(p,'a_hkl0',default, ...
    @(x) validateattributes((x),{'numeric'},{'scalar'}));
%% b_hkl0
default = 0.5;
addParameter(p,'b_hkl0',default, ...
    @(x) validateattributes((x),{'numeric'},{'scalar'}));
%% t_hkl0
default = mean(tof);
addParameter(p,'t_hkl0',default, ...
    @(x) validateattributes((x),{'numeric'},{'scalar'}));
%% sigma0
default =   5*abs(tof(2)-tof(1)); % width
addParameter(p,'sigma0',default, ...
    @(x) validateattributes((x),{'numeric'},{'scalar','positive'}));
%% tau0
default = 5*abs(tof(2)-tof(1)); % assymetry ;
addParameter(p,'tau0',default, ...
    @(x) validateattributes((x),{'numeric'},{'scalar','positive'}));
%% starRange
default = [tof(1) tof(30)];
addParameter(p,'startRange',default, ...
    @(x) validateattributes((x),{'numeric'},{'vector','numel',2,'increasing','>=',tof(1),'<=',tof(end)}));
%% endRange
default = [tof(end-30) tof(end)];
addParameter(p,'endRange',default, ...
    @(x) validateattributes((x),{'numeric'},{'vector','numel',2,'increasing','>=',tof(1),'<=',tof(end)}));
switch lower(p.Results.method)
    case 'attenuation'
        
    case 'attenuationalt'
        %% alpha0
        default = 1e3;
        addParameter(p,'alpha0',default, ...
            @(x) validateattributes((x),{'numeric'},{'scalar'}));
        %% beta0
        default = 1e3;
        addParameter(p,'beta0',default, ...
            @(x) validateattributes((x),{'numeric'},{'scalar'}));
    case '5param'
        default = tof(1) + (tof(end)-tof(1))*[0.33 0.66];
        addParameter(p,'range',default, ...
            @(x) validateattributes((x),{'numeric'},{'vector','numel',2,'increasing','>=',tof(1),'<=',tof(end)}));
    case {'gp','gpcc','gpcc2'}
        %% opts.sig_f  = 1;            %Squared-Exponential Kernel Hyperparameter, output variance
        default = 1;
        addParameter(p,'sig_f',default, ...
            @(x) validateattributes((x),{'numeric'},{'scalar','positive'}));
        %% opts.l      = 1e-4;         %Squared-Exponential Kernel Hyperparameter, lengthscale
        default = 2e-1;
        addParameter(p,'l',default, ...
            @(x) validateattributes((x),{'numeric'},{'scalar','positive'}));
        %% opts.ns     = 3000;         %Number of MC samples used to estimate bragg-edge location and variance.
        default = 300;
        addParameter(p,'ns',default, ...
            @(x) validateattributes((x),{'numeric'},{'scalar','positive','integer'}));
        %% opts.n      = 2500;         %Number of points to sample the Bragg-Edge function.
        default = 3000;
        addParameter(p,'n',default, ...
            @(x) validateattributes((x),{'numeric'},{'scalar','positive','integer'}));
        %% opts.GPscheme   = 'hilbertspace';   %
        default = 'hilbertspace';
        addParameter(p,'GPscheme',default, ...
            @(x) validateattributes((x),{'numeric'},{'scalar','positive','integer'}));
        %% opts.covfunc     = 'M52';
        default = 'M52';
        addParameter(p,'covfunc',default, ...
            @(x) any(validatestring(x,{'M52','M32','SE'})));
        %% opts.optimiseHP     = '';
        default = 'all';
        addParameter(p,'optimiseHP',default, ...
            @(x) any(validatestring(x,{'M52','M32','SE'})));
    case 'crosscorr'
        %% range
        default = tof(1) + (tof(end)-tof(1))*[0.25 0.75];
        addParameter(p,'range',default, ...
            @(x) validateattributes((x),{'numeric'},{'vector','numel',2,'increasing','>=',tof(1),'<=',tof(end)}));
        %% opts.order      = 3;            %polynomial order of Savitzky-Golay
        default = 3;
        addParameter(p,'order',default, ...
            @(x) validateattributes((x),{'numeric'},{'scalar','positive','integer'}));
        %% opts.frame      = 9;            %frame-length for Savitzky-Golay (must be Odd)
        default = 9;
        addParameter(p,'frame',default, ...
            @(x) validateattributes((x),{'numeric'},{'scalar','positive','odd'}));
        %% opts.peakWindow = 20;           %window for peak fitting
        default = 20;
        addParameter(p,'peakWindow',default, ...
            @(x) validateattributes((x),{'numeric'},{'scalar','positive','integer'}));
        %% Peak
        %% opts.y0  
        default = -1;
        addParameter(p,'y0',default, ...
            @(x) validateattributes((x),{'numeric'},{'scalar'}));
        %% opts.xc
        default = 0;
        addParameter(p,'xc',default, ...
            @(x) validateattributes((x),{'numeric'},{'scalar'}));
        %% opts.A
        default = exp(4.5);
        addParameter(p,'A',default, ...
            @(x) validateattributes((x),{'numeric'},{'scalar'}));
        %% opts.Wg
        default = exp(4.3);
        addParameter(p,'Wg',default, ...
            @(x) validateattributes((x),{'numeric'},{'scalar'}));
        %% opts.Wl
        default = exp(4.38);
        addParameter(p,'Wl',default, ...
            @(x) validateattributes((x),{'numeric'},{'scalar'}));
        %% opts.Mu
        default = exp(-3.5);
        addParameter(p,'Mu',default, ...
            @(x) validateattributes((x),{'numeric'},{'scalar'}));

end
parse(p,method,varargin{:});
opts = p.Results;

end