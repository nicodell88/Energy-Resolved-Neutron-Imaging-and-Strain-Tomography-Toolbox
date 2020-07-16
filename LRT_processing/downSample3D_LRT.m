function  [Tr,wl,entry,exit,L,nhat,yInds,nSegs]  = downSample3D_LRT(OB_pro,Proj_pro,sample,mask,rODn,Rno,opts)
%downSample2DLRT averages time of flight data over columns if pixels to
%produce a set of data ready which can subsequently have edges fit to it
%and used for tomography.
%[Tr,wl,entry,exit,L,nhat,yInds,nSegs] = downSample2D_LRT(OB,Proj,rBOo,rODn,Rno,opts)
%
%   Inputs:
%       - OB_pro - is a structure containing:
%           OB_pro.stack :  is the normalised time of flight data for the open beam.
%           OB_pro.lambda:  is the wavelengths associated with each
%                           time-of-flight
%       - Proj_pro - is a structure containing:
%           Proj_pro.stack   :  is the normalised time of flight data for a projection.
%           Proj_pro.lambda  :  is the wavelengths associated with each
%                               time-of-flight.
%       - sample is a structure containing:
%           N: normals to each face, N(i,:) is the normal to the ith face
%           C: centers of each face, C(i,:) is the center of the ith face
%           F: Indices of vertices corresponding to each face, F(i,:) are the
%               indices of the vertices corresponding to face i
%           V: vertices of the shape, V(k,:) = [x,y,z] of the kth vertices
%       - mask is a 512-by-512 matrix of where NaN's indicate pixels
%           outside the sample.
%       - rODn is a 3-by-1 vector defining the sample origin (O) with
%           respect to the centre of the detector (D) in beam coordinates
%           {n}.
%       - Rno is a rotation matrix defining the rotation between beam
%           coordinates {n} and sample coordiantes {o}, such that
%           Rno.'*rABn would rotate the vector rABn from beam coordinates
%           {n} to sample coordinates {o}, producing rABo.
%       - opts is a structure containing:
%           opts.nPix : specifies to use nPix-by-nPix pixels for
%               downsampling.
%           opts.nanThresh indicated what fraction of pixels must lie
%               within the mask for the measurement to be counted.
%
% See also GenerateLRT3D, find_intersects_3D, GenerateLRT2D.


% Copyright (C) 2020 The University of Newcastle, Australia
% Authors:
%   Nicholas O'Dell <Nicholas.Odell@newcastle.edu.au>
% Last modified: 16/07/2020
% This program is licensed under GNU GPLv3, see LICENSE for more details.


%% Downsample
nPixRow = size(Proj_pro.stack,1);
nPixCol = size(Proj_pro.stack,2);


nres = size(Proj_pro.stack,1);
nwl  = numel(Proj_pro.lambda);
PI = NaN(nres^2,1);
PJ = NaN(nres^2,1);
PAve = NaN(nres^2,nwl);
OBAve = NaN(nres^2,nwl);

iter=0;

try
    h = waitbar(0,'Downsampling');
    Hfig = figure(1);
    totalAtten = mean(Proj_pro.stack,3);
    clf
    H =  pcolor(totalAtten);
    shading flat
    
    
    for i =1:(floor(nPixRow/opts.nPix))
        for j     =1:(floor(nPixCol/opts.nPix))
            X =  ((i-1)*(floor(nPixCol/opts.nPix))+j)/...
                (floor(nPixRow/opts.nPix)*floor(nPixCol/opts.nPix));
            waitbar(X,h);
            
            i_inds = (i-1)*opts.nPix +(1:(opts.nPix));% - round(opts.nPix/2);
            j_inds = (j-1)*opts.nPix +(1:(opts.nPix));% - round(opts.nPix/2);
            
            idx = i_inds>0 & i_inds<=nPixRow;
            i_inds = i_inds(idx);
            
            idx = j_inds>0 & j_inds<=nPixCol;
            j_inds = j_inds(idx);
            
            
            n_pix = length(i_inds)*length(j_inds);
            
            
            
            if sum(~isnan(mask(i_inds,j_inds)),'all') ...
                    /numel(mask(i_inds,j_inds)) < (opts.nanThresh)
                continue
            else
                iter=iter+1;
                [I,J] = meshgrid(i_inds,j_inds);
                I = I(:);
                J = J(:);
                
                %avoiding for loop...
                a = 1:nwl;
                b = repelem(a',numel(I));
                newJ = repmat(J,nwl,1);
                newI = repmat(I,nwl,1);
                
                inds = sub2ind(size(Proj_pro.stack),newI,newJ,b);
                indsplot = sub2ind(size(totalAtten),newI,newJ);
                
                
                totalAtten(indsplot) = nan;
                set(H,'CData',totalAtten);
                
                
                Proj_sec = Proj_pro.stack(inds);
                OB_sec = OB_pro.stack(inds);
                Proj_sec = permute(Proj_sec,[2,1,3]);
                OB_sec = permute(OB_sec,[2,1,3]);
                Proj_sec = reshape(Proj_sec,n_pix,nwl);
                OB_sec = reshape(OB_sec,n_pix,nwl);
                
                
                PAve(iter,:) = squeeze(nanmean(Proj_sec,1));
                OBAve(iter,:) = squeeze(nanmean(OB_sec,1));
                PI(iter) = mean(I(~isnan(Proj_sec(:,1))));
                PJ(iter) = mean(J(~isnan(Proj_sec(:,1))));
                

            end
        end
    end
catch me
    close(h)
    close(Hfig)
    rethrow(me)
end
close(h)
close(Hfig)
%% Delete stuff
OB_cell = OBAve(~isnan(PI),:);
Proj_cell = PAve(~isnan(PI),:);
I_cell =  PI(~isnan(PI));
J_cell = PJ(~isnan(PI));
% Calculate the transmission while we're here....
Tr = [Proj_cell]./[OB_cell];

%% Geometric Properties
rPDn = pix2vec(I_cell,J_cell);

%map detector to sample
[~,n] = size(rPDn);
nhatt = repmat([1;0;0],1,n);
nhat_o  = Rno.'*nhatt;
rPOo = Rno.'*(rPDn - rODn);

lines = [rPOo;nhat_o];

[entry,exit,nhat,L,yInds,nSegs,rejects] = find_intersects_3D(sample,lines);


wl = Proj_pro.lambda;
end

