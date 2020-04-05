function [StrainImage,SigmaImage,opts] = makeStrainImage(OB,Proj,opts,d0Tr)
%MAKESTRAINIMAGE Generates a strain-image from a single projection
%   [StrainImage,SigmaImage,opts] = makeStrainImage(OB,Proj,opts)
%   Inputs:
%       - OB is a structure containing data from the open beam projection.
%           OB.im_stack     :   a 3D array of time of flight transmission
%                               data. M-by-M-by-N where N is the number of
%                               wavelengths recorded.
%           OB.ntrig        :   a 1-by-N array of triggers recorded per
%                               wavelength.
%           OB.tof          :   a N-by-1 array of wavelengths.
%       - Proj is a structure containing the data from a projection, It has
%       the same fields as OB.
%       - opts is a structure containing
%           opts.supplyMask :   a char array defining the method used to
%                               supply the pixel mask defining which rays
%                               passed through the material. ('gui','auto',
%                               'supply'). AUTO attempts to automatically
%                               generate the mask using edge height as a
%                               metric. GUI allows the user to select the
%                               verticies of the "shadow" cast on the
%                               detector by the beam.
%           opts.rangeLeft  :   a 2 element vector containing the start and
%                               end wavelengths for determining the edge
%                               height on the left side of the edge.
%           opts.rangeRight :   a 2 element vector containing the start and
%                               end wavelengths for determining the edge
%                               height on the right side of the edge.
%           opts.d0         :   time-of-flight corresponding to the
%                               unstrained lattice parameter.
%           opts.nPix       :   Size of macro-pixels for averaging over.
%           opts.Thresh     :   Maximum fraction of pixels allowed within
%                               macro-pixel but not within mask. e.g., 0.05
%           opts.maskThresh :   Threshold value for determining mask...
%                               TODO: Fix description.
%           opts.sigma_d0   :   Standard deviation of the supplied d0;
%           opts.BraggOpts  :   A structure containing all of the options
%
%   Outputs:
%       - StrainImage is a 2D array containing the strain measurements,to
%       be plotted as an image.
%       - SigmaImage is a 2D array, containing the estimated standard
%       deviation of the
%       - opts is the modified options structure, with the same form as the
%       input
%
% Copyright (C) 2020 The University of Newcastle, Australia
% Authors:
%   Nicholas O'Dell <Nicholas.Odell@newcastle.edu.au>
% Last modified: 06/03/2020
% This program is licensed under GNU GPLv3, see LICENSE for more details.

%% process inputs
if ~isfield(opts,'supplyMask')
    opts.supplyMask = 'auto';
end

assert(isfield(opts,'nPix'),...
    'opts.npix, defines the size of the macropixels and must be defined.')

assert(isfield(opts,'d0'),...
    'opts.d0, defines the unstrained lattice parameter and must be defined.')

%% Normalise OB and Projection=
%OB
ntrigs_rep_OB = reshape(OB.ntrigs,1,1,numel(OB.ntrigs)).*ones(size(OB.im_stack,1),size(OB.im_stack,2),numel(OB.ntrigs));
OB.im_stack = OB.im_stack./ntrigs_rep_OB;
%Proj
ntrigs_rep_Proj = reshape(Proj.ntrigs,1,1,numel(Proj.ntrigs)).*ones(size(Proj.im_stack,1),size(Proj.im_stack,2),numel(Proj.ntrigs));
Proj.im_stack = Proj.im_stack./ntrigs_rep_Proj;

Tr = Proj.im_stack ./ OB.im_stack;
%% Create/Obtain mask
switch lower(opts.supplyMask)
    case 'auto'
        [~,opts.startIdx] = min((OB.tof(:).'-opts.rangeLeft(:)).^2,[],2);
        [~,opts.endIdx] = min((OB.tof(:).'-opts.rangeRight(:)).^2,[],2);
        
        edge = log(mean(Tr(:,:,opts.endIdx(1):opts.endIdx(2)),3))-log(mean(Tr(:,:,opts.startIdx(1):opts.startIdx(2)),3));
        edge(isinf(edge)) = 0;
        opts.mask = edge>opts.maskThresh;
    case 'supply'
        assert(isfield(opt,'mask'),'Mask must be supplied if opts.supplyMask = ''supply''')
    case 'gui'
        %% Obtain mask from user
        [~,opts.startIdx] = min((OB.tof(:).'-opts.rangeLeft(:)).^2,[],2);
        [~,opts.endIdx] = min((OB.tof(:).'-opts.rangeRight(:)).^2,[],2);
        
        edge = log(nanmean(Tr(:,:,opts.endIdx(1):opts.endIdx(2)),3))-log(nanmean(Tr(:,:,opts.startIdx(1):opts.startIdx(2)),3));
        edge(isinf(edge)) = 0;
        
        figure(1)
        clf
        pcolor(edge)
        shading flat
        daspect([1 1 1])
        msg = sprintf('Select the verticies of the sample in a\nclockwise direction the image.\nPress ENTER when finished');
        title(msg)
        [xv,yv]=ginput();
        %% Calculate Mask
        idx = 1:numel(Proj.im_stack(:,:,1));
        [yq,xq] = ind2sub(size(Proj.im_stack(:,:,1)),idx);
        idx = inpolygon(xq,yq,[xv;xv(1)],[yv;yv(1)]);
        
        opts.mask = zeros(size(Proj.im_stack(:,:,1)));
        opts.mask(idx) = 1;
        opts.mask = logical(opts.mask);
        
        %% display Mask
        imagesc(opts.mask)
        daspect([1 1 1])
        title('Mask')
        xlabel('X - [pixels]')
        ylabel('Y - [pixels]')
    otherwise
        error('opts.supplyMask must be ''auto'', ''supply'' or ''gui''.')
end

assert(isfield(opts,'mask'),'User supplied mask must be specified is ''''opts.supplyMask'''' is true.')
assert(islogical(opts.mask),'opts.mask must be a MxN array of logicals where M and N are the number of pixels vertically and horizontally respectively.')


tmp = opts.mask;
opts.mask = nan(size(tmp));
opts.mask(tmp) = 1;
%% Downsample
nPixRow = size(Proj.im_stack,1);
nPixCol = size(Proj.im_stack,2);

indicator_macro = zeros(length(1:(floor(nPixRow/opts.nPix))),length(1:(floor(nPixCol/opts.nPix))));

nwl = length(Proj.tof);

Proj_masked = Proj.im_stack .*opts.mask;
OB_masked   = OB.im_stack   .*opts.mask;

PAve = [];
OBAve = [];
PI = [];

for j = 1:(floor(nPixCol/opts.nPix))        %Order of for loops is important due to mixed indexing
    for i = 1:(floor(nPixRow/opts.nPix))
        
        i_inds = (i-1)*opts.nPix +(1:(opts.nPix));
        j_inds = (j-1)*opts.nPix +(1:(opts.nPix));
        
        n_pix = length(i_inds)*length(j_inds);
        
        indicator_macro(i,j)    = sum(~isnan(opts.mask(i_inds,j_inds)),'all') ...
            /numel(opts.mask(i_inds,j_inds)) > (1-opts.Thresh);
        
        if ~indicator_macro(i,j)
            continue
        else
            [I,J] = meshgrid(i_inds,j_inds);
            I = I(:);
            J = J(:);
            
            %avoiding for loop...
            a = 1:nwl;
            b = repelem(a',numel(I));
            newJ = repmat(J,nwl,1);
            newI = repmat(I,nwl,1);
            inds = sub2ind(size(Proj_masked),newJ,newI,b);
            Proj_sec = Proj_masked(inds);
            OB_sec = OB_masked(inds);
            Proj_sec = permute(Proj_sec,[2,1,3]);
            OB_sec = permute(OB_sec,[2,1,3]);
            Proj_sec = reshape(Proj_sec,n_pix,nwl);
            OB_sec = reshape(OB_sec,n_pix,nwl);
            
            
            PAve    = [PAve;squeeze(nanmean(Proj_sec,1))];
            OBAve   = [OBAve;squeeze(nanmean(OB_sec,1))];
            tmp = mean(I(~isnan(Proj_sec(:,1))));
            
            if isnan(tmp)
                indicator_macro(i,j) = 0;
            end
            PI = [PI;tmp];
            
        end
        
    end
end
OB_cell = {OBAve(~isnan(PI),:)};
Proj_cell = {PAve(~isnan(PI),:)};

Tr_cell = {[Proj_cell{1}]./[OB_cell{1}]};

%% Fit Bragg Edges
if exist('d0Tr','var')
    [d_cell,std_cell,~,~,opts.BraggOpts] = fitEdges(Tr_cell,Proj.tof,opts.BraggOpts,d0Tr);
else
    [d_cell,std_cell,~,~,opts.BraggOpts] = fitEdges(Tr_cell,Proj.tof,opts.BraggOpts);
end
%% Shuffle Data, produce plots
idx = find(logical(indicator_macro));
StrainImage = nan(size(indicator_macro));
SigmaImage  = nan(size(indicator_macro));

if strcmpi(opts.BraggOpts.method,'crosscorr')
    StrainImage(idx) = (d_cell{1})/opts.d0;
else
    StrainImage(idx) = (d_cell{1}-opts.d0)/opts.d0;
end

% SigmaImage(idx)  = std_cell{1}/opts.d0; %cheating

if ~isfield(opts,'sigma_d0')
    warning('A standard deviation for d0 has not been supplied and the minimum of the confidence intervals from the fitting procedure has been used.')
    sigmad0 = min(std_cell{1});
else
    sigmad0 = opts.sigma_d0;
end

SigmaImage(idx) = sqrt(...
    (1/opts.d0)^2*std_cell{1}.^2 + (-d_cell{1}*(opts.d0^(-2))).^2 * sigmad0^2 ...
    );

end