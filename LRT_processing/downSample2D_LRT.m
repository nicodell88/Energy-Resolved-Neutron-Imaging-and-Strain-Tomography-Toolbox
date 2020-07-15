function [Tr,wl,entry,exit,L,nhat,yInds,nSegs] = downSample2D_LRT(OB_pro,Proj_pro,rBOo,mask,rODn,Rno,opts)
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
%       - rBOo is a 2-by-N matrix, where each column is an x,y vertex
%           of the sample boundary. Closed polygons, annular shapes and
%           multi-body samples can be specified. E.g., a closed polygon, in 
%           this case a square, can be specified as:
%               rBOo = [0,0,1,1,0;0,1,1,0,0]; 
%           Annular or multi body shape boundaries are expressed by placing
%           a column of NaNs between the rows defining each shape. E.g., a
%           square within a square;
%               rBOo = [0,0,1,1,0,NaN,0.1,0.1,0.9,0.9,0.1;
%                       0,1,1,0,0,NaN,0.1,0.9,0.9,0.1,0.1];
%       - mask is a 512-by-1 vector of logicals indicating which columns
%           are associated with pixels that lie behind the sample.
%       - rODn is a 2-by-1 vector defining the sample origin (O) with
%           respect to the centre of the detector (D) in beam coordinates
%           {n}.
%       - Rno is a rotation matrix defining the rotation between beam
%           coordinates {n} and sample coordiantes {o}, such that
%           Rno.'*rABn would rotate the vector rABn from beam coordinates
%           {n} to sample coordinates {o}, producing rABo.
%       - opts is a structure containing:
%           opts.nCol is the number of columns to average over.
%           opts.rowRange is a 2 element vector defining the start and end
%           row to average over.
%
% See also GenerateLRT2D, find_intersects_2D, GenerateLRT3D.

% Copyright (C) 2020 The University of Newcastle, Australia
% Authors:
%   Nicholas O'Dell <Nicholas.Odell@newcastle.edu.au>
% Last modified: 15/07/2020
% This program is licensed under GNU GPLv3, see LICENSE for more details.

%% Downsample
nres = size(Proj_pro.stack,1);
nwl  = numel(Proj_pro.lambda);
PI = NaN(nres,1);
PAve = NaN(nres,nwl);
OBAve = NaN(nres,nwl);
indicator_macro = zeros(nres,1);

try
    h = waitbar(0,'Downsampling');
    Hfig = figure(1);
    totalAtten = mean(Proj_pro.stack,3);
    clf
    H =  pcolor(totalAtten);
    shading flat
    
    for i =1:(floor(nres/opts.nCol))
        iter=i;
        X = i/(floor(nres/opts.nCol));
        waitbar(X,h);
        j_inds = (i)*opts.nCol + [0:(opts.nCol-1)];
        %         j_inds = (i)*opts.nCol + [-2:+];
        i_inds = opts.rowRange(1):opts.rowRange(2);
        
        idx = i_inds>0 & i_inds<=nres;
        i_inds = i_inds(idx);
        
        idx = j_inds>0 & j_inds<=nres;
        j_inds = j_inds(idx);
        
        %         if(~all(opts.AlignResult.mask(j_inds,k)==1))
        %
        %         end
        n_pix = length(i_inds)*length(j_inds);
        
        indicator_macro(i)    = all(mask(j_inds)==1);
        
        if ~indicator_macro(i)
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
            %                         inds = sub2ind(size(Proj_pro.stack),newJ,newI,b);
            inds = sub2ind(size(Proj_pro.stack),newI,newJ,b);
            indsplot = sub2ind(size(totalAtten),newI,newJ);
            %             indsplot = sub2ind(size(totalAtten),newJ,newI);
            
            totalAtten(indsplot) = nan;
            set(H,'CData',totalAtten);
            
            
            Proj_sec = Proj_pro.stack(inds);
            OB_sec = OB_pro.stack(inds);
            Proj_sec = permute(Proj_sec,[2,1,3]);
            OB_sec = permute(OB_sec,[2,1,3]);
            Proj_sec = reshape(Proj_sec,n_pix,nwl);
            OB_sec = reshape(OB_sec,n_pix,nwl);
            
            
            %             PAve    = [PAve;squeeze(nanmean(Proj_sec,1))];
            %             OBAve   = [OBAve;squeeze(nanmean(OB_sec,1))];
            PAve(iter,:)    = squeeze(nanmean(Proj_sec,1));
            OBAve(iter,:)   = squeeze(nanmean(OB_sec,1));
            
            tmp = mean(J(~isnan(Proj_sec(:,1))));
            
            if isnan(tmp)
                indicator_macro(i) = 0;
            end
            PI(iter) = tmp;
        end
    end
catch me
    close(h)
    close(Hfig)
    rethrow(me)
end
close(h)
close(Hfig)

% append to structure
% Delete the extra NaNs and store...

OB_cell = OBAve(~isnan(PI),:);
Proj_cell = PAve(~isnan(PI),:);
I_cell =  PI(~isnan(PI));
% Calculate the transmission while we're here....
Tr = [Proj_cell]./[OB_cell];

%% Geometric properties

rPDn = [zeros(1,length(I_cell));
    55e-6*(I_cell(:).'-513/2)];%55e-6*(v-513/2)

%map detector to sample
[~,n] = size(rPDn);
nhatt = repmat([1;0],1,n);
nhat_o  = Rno.'*nhatt;
rPOo = Rno.'*(rPDn - rODn);


[entry,exit,nhat,L,yInds,nSegs,Rejects] = find_intersects_2D(rBOo,rPOo,nhat_o);

wl = Proj_pro.lambda;
end