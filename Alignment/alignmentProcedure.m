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
% Last modified: 08/05/2020
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

%% read STL
%remember to convert dimensions
%% create FE model
%remember to convert dimensions
%% Downsample nodes and shadow

%% Set up optimisation options

%% run alignment

%% convert results to absolute sample position

%% produce mask

rODn = 0;
Rno = 0;
mask = 0;
X = 0;



end

