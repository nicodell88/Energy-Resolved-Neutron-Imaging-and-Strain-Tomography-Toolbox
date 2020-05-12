function [rODn,Rno,mask,X] = alignmentProcedure(stlFile,stlUnit,rHSs_all,Rsh_all,Shadow,varargin)
%[rODn,Rno,mask,X] = alignmentProcedure(stlFile,stlUnit,rHSs_all,Rsh_all,varargin)
% Inputs:
%   - stlFile is a char-array or string containing a relative or absolute
%       filepath to the stl file of the sample.
%   - stlUnits is a char-array or string which indicates the unit of the
%       CAD model, must be {'mm','m','inch'}
%   -rHSs_all is a 3xN matrix
%   - opts is a structure containing
%       opts.rangeLeft  :   a 2 element vector containing the window for
%                           averaging the edge height before the edge
%       opts.rangeRight :   a 2 element vector containing the window for
%                           averaging the edge height after the edge
%
% Outputs:
%   - edges is a mPix-by-nPix-by-nProj matrix where each page is contains the
%     "edge heights" for each projection
%
% Copyright (C) 2020 The University of Newcastle, Australia
% Authors:
%   Nicholas O'Dell <Nicholas.Odell@newcastle.edu.au>
% Last modified: 12/05/2020
% This program is licensed under GNU GPLv3, see LICENSE for more details.

p = inputParser;

%% STL
addRequired(p,'stlFile',@(x) assert(exist(x,'file') == 2,'''stlFile'' must exist'));
addRequired(p,'stlUnit',@(x) any(validatestring(x,{'mm','m','inch'})));
%% Known Transform
validate_rHSs_all = @(x) size(x,1)==3 && isnumeric(x);
addRequired(p,'rHSs_all',validate_rHSs_all)
validate_Rsh_all = @(x) size(x,1)==3 && size(x,2)==3 && isnumeric(x);
addRequired(p,'Rsh_all',validate_Rsh_all)
%% Shadow
addRequired(p,'Shadow',@(x) validateattributes(x,{'numeric','logical'},{'ndims',3,'size',[512,512,nan]}))
%% Optional Params
addParameter(p,'sigma'          ,3e-4 ,@(x) validateattributes(x,{'numeric'},{'positive'}));
addParameter(p,'plotGrid'       ,[3,4]  ,@(x) (validateattributes(x,{'numeric'},{'vector','positive','integer','numel',2})));
addParameter(p,'Nnodes'         ,1.5e3  ,@(x) (validateattributes(x,{'numeric'},{'scalar','integer','>',1e3})));
addParameter(p,'Npix'           ,1.5e3  ,@(x) (validateattributes(x,{'numeric'},{'scalar','integer','>',1e3})));
addParameter(p,'thresh'         ,0      ,@(x) validateattributes(x,{'numeric'},{'scalar'}) );       
%% Parse inputs
parse(p,stlFile,stlUnit,rHSs_all,Rsh_all,Shadow,varargin{:})
%% Check dimensions match
assert(isequal(size(p.Results.Rsh_all,3),size(p.Results.rHSs_all,2), size(p.Results.Shadow,3)),'Expected rHSs_all to be 3xN, Rhs_all to be 3x3xN, and Shadow to be 512x512xN');
%% Conversion factor
switch lower(p.Results.stlUnit)
    case 'm'
        convert2m = 1;
    case 'mm'
        convert2m = 1e-3;
    case 'inch'
        convert2m = 25.4e-3;
end
%% read STL
[sample.F,sample.V] = stlread(p.Results.stlFile);
sample.V = sample.V*convert2m;
%remember to convert dimensions
%% create FE model
%remember to convert dimensions
model = createpde(1);
importGeometry(model,p.Results.stlFile);

hMax = 10;

while true
    mesh = generateMesh(model,'Hmax',hMax);
    nNodes = length(mesh.Nodes);
nNodes_des = p.Results.Nnodes;

if nNodes >= nNodes_des
   break 
end

hMax = hMax*0.95;

end

rVOo = mesh.Nodes*convert2m;

%% Downsample shadow
np = size(p.Results.Shadow,3);
rPDn = cell(np,1);
rPDnDownsampled = cell(np,1);
beta = cell(np,1);
for i = 1:np
    patch = p.Results.Shadow(:,:,i);
    patch(isinf(patch)) = 0;
    
    [Nrow,Ncol] = size(patch);
    centres = pix2vec(1:Nrow,1:Ncol);
    yCentres = centres(2,:);
    zCentres = centres(3,:);
    [Y,Z] = meshgrid(yCentres.',zCentres.');
    
    idx   = find(patch>p.Results.thresh);
    [u,v] = ind2sub(size(patch),idx);
    rPDn{i} = pix2vec(u,v);
    
    IDX = find(patch>p.Results.thresh);
    idx = randsample(IDX,p.Results.Npix,false);
    rPDnDownsampled{i} = [zeros(length(idx),1).';Y(idx).';Z(idx).'];
    
    beta{i} = ones(length(rPDnDownsampled{i}),1);
end
%% Set up optimisation options
rPDd = rPDnDownsampled(1:end);
Theta_ns    = [0,0,0];
Theta_ho    = [0,0,0];
rSDn        = [0,0];
rOHh        = [0,0,0];

x0 = [...
    Theta_ns(:);
    Theta_ho(:);
    rSDn(:);
    rOHh(:)];

opts.sigma = p.Results.sigma;
opts.plotGrid = p.Results.plotGrid;

plotFcn = @(x,state,vars) plotProj(x,state,sample,p.Results.Shadow,rHSs_all,Rsh_all,opts,vars);

options = optimoptions(@fmincon,...
    ...'algorithm','sqp',...
    'algorithm','trust-region-reflective',...
    'specifyObjectiveGradient',true,...
    'checkGradients',false,...
    'display','iter',...
    'FunctionTolerance',1e-8,...
    'SubProblemAlgorithm','factorization',...
    'HessianFcn','objective',...
    'OptimalityTolerance',2e4,...
    'StepTolerance',1e-8,...
    'MaxIterations',20,...
    'OutputFcn',plotFcn,...
    'PlotFcn',{@optimplotx,@optimplotfval});

% beta = ones();
costFun = @(x)AlignmentCostFunction([x],rPDd,beta,rVOo,rHSs_all,Rsh_all,opts);

%% Testing
state.iteration = 1;
vars = 'iter';
plotFcn(zeros(11,1),state,vars)
%% run alignment
Xopt = fmincon(costFun,x0,[],[],[],[],[],[],[],options);
%% convert results to absolute sample position

%% produce mask

rODn = 0;
Rno = 0;
mask = 0;
X = 0;



end

