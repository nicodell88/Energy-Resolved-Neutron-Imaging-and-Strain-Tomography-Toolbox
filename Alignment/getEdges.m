function edges = getEdges(projPath,opts)
%getEdges(projPath,opts) calculates the "pseudo edge heights" for use in
%   sample alignment.
% Inputs:
%   - projPath is a relative filepath to the folder containing the
%   projections "im_stack". XXX_proj_n.mat
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
% Last modified: 07/04/2020
% This program is licensed under GNU GPLv3, see LICENSE for more details.

% Grab data files containing "proj"
% pathStr = fullfile(projPath,'*proj*.mat');
% dataFiles = dir(pathStr);
% assert(length(dataFiles)>=1,'Expected more than zero projections to exist!');
OB = load(opts.OB,'tof','im_stack','ntrigs');
ntrigs_rep_OB = reshape(OB.ntrigs,1,1,numel(OB.ntrigs)).*ones(size(OB.im_stack,1),size(OB.im_stack,2),numel(OB.ntrigs));
OB.im_stack = OB.im_stack./ntrigs_rep_OB;

fname = sprintf(opts.fmt,opts.proj_idx(1));
dataFiles = dir(fullfile(projPath,fname));

% Determine data size
matObj = matfile(fullfile(dataFiles(1).folder,dataFiles(1).name));

nPixRow = size(matObj,'im_stack',1);
nPixCol = size(matObj,'im_stack',2);
nProj     = length(opts.proj_idx);

% Initialise edges
edges = nan(nPixRow,nPixCol,nProj);


delete(findall(0,'tag','TMWWaitbar'));
wh      = waitbar(0,'calculating edge height', ...
    'Name', 'Bragg Edge Progress Bar', ...
    'CreateCancelBtn', 'setappdata(gcbf,''cancelling'',1)');
for i = 1:nProj
    waitbar(i/nProj,wh); % Update waitbar
    if getappdata(wh,'cancelling')
        warning('User cancelled operation');
        break
    end
    
    str = sprintf(opts.fmt,opts.proj_idx(i));
    
    Proj = load(fullfile(projPath,str),'im_stack','tof','ntrigs');
    ntrigs_rep_Proj = reshape(Proj.ntrigs,1,1,numel(Proj.ntrigs)).*ones(size(Proj.im_stack,1),size(Proj.im_stack,2),numel(Proj.ntrigs));
    Proj.im_stack = Proj.im_stack./ntrigs_rep_Proj;
    
    [~,leftIdx]    = min((Proj.tof(:).'-opts.rangeLeft(:)).^2,[],2);
    [~,rightIdx]   = min((Proj.tof(:).'-opts.rangeRight(:)).^2,[],2);
    
    Iright  = mean(Proj.im_stack(:,:,rightIdx(1):rightIdx(2)),3);
    ILeft   = mean(Proj.im_stack(:,:,leftIdx(1):leftIdx(2)),3);
    I0left  = mean(OB.im_stack(:,:,leftIdx(1):leftIdx(2)),3);
    I0right = mean(OB.im_stack(:,:,rightIdx(1):rightIdx(2)),3);
    
    edges(:,:,i) = log(Iright) -log(I0right) + log(I0left) - log(ILeft);
    
    
    %     before_edge = mean(Proj.im_stack(:,:,leftIdx(1) :leftIdx(2) ),3);
    %     after_edge  = mean(Proj.im_stack(:,:,rightIdx(1):rightIdx(2)),3);
    
    %    edge = log(mean(Tr(:,:,opts.endIdx(1):opts.endIdx(2)),3))-log(mean(Tr(:,:,opts.startIdx(1):opts.startIdx(2)),3));
    %     edges(:,:,i) = log(after_edge) - log(before_edge);
    
end
delete(wh);
end