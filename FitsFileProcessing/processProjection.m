function [data] = processProjection(fileNames,pathName,opts)
% [data] = processSingleFitFile(fileNames,opts.path)
%   Processes the fits files and returns a 3d matrix containing the
%   intensities
%   Operates on a single projection.
%
% Inputs:
%   - fileNames: cell array containing the names of the fits files for this
% one projection.
%	- pathName: path to the directory containing the names of the fits files to be processed.
%
% Outputs:
%   - data is a structure containing
%   	- im_stack		:	a 3d hypermatrix where each value corresponds to the detector counts
%   						of a pixel.
%   						The first two dimensions of the matrix are spatial and correspond to the physical pixel positions.
%   						The third dimension is temporal (i.e. wavelengths).
%
%   	- tof 			:	The time of flight corresponding to the temporal direction in
%   image_stack. i.e the TOF for each (:,:,i).
%
%   	- total_counts	: the total counts per image frame
%
% Copyright (C) 2020 The University of Newcastle, Australia
% Authors:
%   Nicholas O'Dell <Nicholas.Odell@newcastle.edu.au>
% Last modified: 21/04/2020
% This program is licensed under GNU GPLv3, see LICENSE for more details.

nn = length(fileNames);
%% Initialise the structure

% im_stack
% 512 pixels x 512 pixels x length(idx) wavelengths.
% each (:,:,i) is a wavelength.
% each (i,j,~) is the counts at a pixel.
data.im_stack = zeros(512,512,nn);

% tof
% length(idx) wavelengths x 1.
% The time-of-flight for each wavelength (slice in im_stack).
data.tof = zeros(nn,1);

% total_counts
% length(idx) wavelengths x 1.
% total_counts(i) = sum(im_stack(:,:,i)).
data.total_counts = zeros(nn,1);

% ntrigs
% length(idx) wavelengths x 1
% each time a pulse happens --- divide by this.
data.ntrigs = nan(nn,1);


%% Loop over files and process

h = waitbar(0,'Processing projection...','interpreter','tex');

for i = 1:nn
    n = i;
    %Print out waitbar with message
    msg = sprintf('Processing: %s',fileNames{n});
    waitbar(i/nn,h,msg);
    
    FitPathI = dir(fullfile(pathName,fileNames{n}));
    assert(numel(FitPathI)==1,'Expected exactly 1 matching file for %s',fullfile(pathName,fileNames{n}))
    fileI = fullfile(FitPathI.folder,FitPathI.name);
    info = fitsinfo(fileI);
    %     info=fitsinfo([pathname,fileNames{n}]);
    data.tof(i)=info.PrimaryData.Keywords{9,2};
    data.ntrigs(i) = info.PrimaryData.Keywords{12,2};
    data.total_counts(i)=info.PrimaryData.Keywords{11,2};
    %     data.im_stack(:,:,i) = flipud(fliplr(fitsread(fileI)'));
    data.im_stack(:,:,i) = rot90(fitsread(fileI)',2);
    % transposed to fix the fact taht c code is row major
    % the flip up down and left right put things into our right hand coord system
    
end
delete(h);
data.spatial = mean(data.im_stack,3);
data.idx = opts.tof_idx;


end