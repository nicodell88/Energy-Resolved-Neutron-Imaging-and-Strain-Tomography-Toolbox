function plotStrainImage(StrainImage,SigmaImage,opts)
%PLOTSTRAINIMAGE Plots a strain-image from a single projection.
%   plotStrainImage(StrainImage,SigmaImage,opts)
%   Inputs:
%       - StrainImage is a 2D array containing the strain measurements,to
%       be plotted as an image.
%       - SigmaImage is a 2D array, containing the estimated standard
%       deviation of the
%       - opts is a structure containing.
%           opts.figNum :   Figure number used for plotting.
%           opts.nRes   :   Side length of macro-pixels.
%   Outputs:
%       - None.
%
% Copyright (C) 2020 The University of Newcastle, Australia
% Authors:
%   Nicholas O'Dell <Nicholas.Odell@newcastle.edu.au>
% Last modified: 06/03/2020
% This program is licensed under GNU GPLv3, see LICENSE for more details.
addpath ./utility_functions/

if ~isfield(opts,'figNum')
    opts.figNum = 1;
end


H = figure(opts.figNum);
clf
H.Position = [50 50 1500 600];

H1 = subplot(1,2,1);
tmp = load('RdYlBu_up.mat');
cmap = tmp.RdYlBu;
H = pcolor(StrainImage);
set(H,'alphaData',~isnan(StrainImage));
shading flat
daspect([1 1 1])
xlabel('X - [pixels]')
ylabel('Y - [pixels]')
title('Strain - $\epsilon$')
colorbar
absval = max(abs(StrainImage(:)));
absval = min(absval,1e-3);

caxis([-absval absval])
colormap(gca,cmap)

H2 = subplot(1,2,2);
H = pcolor(SigmaImage);
set(H,'alphaData',~isnan(StrainImage));
shading flat
daspect([1 1 1])
xlabel('X - [pixels]')
ylabel('Y - [pixels]')
title('Standard Deviation - $\sigma$')
colorbar
caxis([0 min(max(SigmaImage(:)),1e-3)])

linkaxes([H1,H2],'xy');

ChangeTicks(H1)
ChangeTicks(H2)


%% Change Ticks
    function ChangeTicks(HH)
        ticks=get(HH,'YTick');%retrieve current ticks
        ticks=ticks*opts.nRes;%multiply
        ticks = num2cell(ticks);
        ticks=cellfun(@num2str,ticks,'UniformOutput',false);%convert to cellstr
        set(HH,'YTickLabels',ticks)%set new tick labels
        
        ticks=get(HH,'XTick');%retrieve current ticks
        ticks=ticks*opts.nRes;%multiply
        ticks = num2cell(ticks);
        ticks=cellfun(@num2str,ticks,'UniformOutput',false);%convert to cellstr
        set(HH,'XTickLabels',ticks)%set new tick labels
        
    end
end