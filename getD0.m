function [Tr,tof] = getD0(OBfile,D0file)
% getD0 generates generates a single edge by averaging over a user
% specified region, the intention of this function is to easily obtain a d0
% edge.
%[Tr,tof] = getD0(OBfile,D0file)
%Inputs:
%   - OBfile is the absolute or relative path to the open beam .mat file.
%   - D0file is the absolute or relative path to the D0 .mat file.
%
%Outputs: 
%   - Tr is a 1xN vector containg the transmission data for the edge.
%   - tof is the time of flight spectra.
%
% See also processFitsFiles



% Copyright (C) 2020 The University of Newcastle, Australia
% Authors:
%   Nicholas O'Dell <Nicholas.Odell@newcastle.edu.au>
% Last modified: 08/07/2020
% This program is licensed under GNU GPLv3, see LICENSE for more details.


OBunPro = load(OBfile);
d0unPro = load(D0file);

addpath ./utility_functions/
OB = ProcessMat(OBunPro);
Proj = ProcessMat(d0unPro);

Tr_ = Proj.stack ./ OB.stack;

%%
figure(1)
clf
pcolor(nanmean(Tr_,3));
shading flat
daspect([1 1 1])
msg = sprintf('Select the verticies of the sample in a\nclockwise direction the image.\nPress ENTER when finished');
title(msg)
[xv,yv]=ginput();
%% Calculate Mask
idx = 1:numel(Proj.stack(:,:,1));
[yq,xq] = ind2sub(size(Proj.stack(:,:,1)),idx);
idx = inpolygon(xq,yq,[xv;xv(1)],[yv;yv(1)]);

mask = zeros(size(Proj.stack(:,:,1)));
mask(idx) = 1;
mask = logical(mask);
%%
ProjMasked = Proj.stack.*mask;
OBMasked = OB.stack.*mask;
tof = OBunPro.tof;
Tr = squeeze(sum(ProjMasked,[1 2])./sum(OBMasked,[1 2]));

end