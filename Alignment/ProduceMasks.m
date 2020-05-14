function xProjection = ProduceMasks(sample, rPDn, rODn,Rno,opts)
% Copyright (C) 2020 The University of Newcastle, Australia
% Authors:
%   Nicholas O'Dell <Nicholas.Odell@newcastle.edu.au>
% Last modified: 14/05/2020
% This program is licensed under GNU GPLv3, see LICENSE for more details.

% Initialise Output
xProjection = nan(opts.nRowPix,opts.nColPix);

rVOo = sample.V;
rVDn = (rODn + Rno*rVOo.').';

% tmp = reshape(rVDn(sample.F(:).',:),3,[]);
% tmp = [tmp;nan(1,length(tmp))];
% poly_rVDn = reshape(tmp,[],3);
% yPix = poly_rVDn(:,2);
% zPix = poly_rVDn(:,3);
F = sample.F.';
tmp = rVDn(F(:),:);

yPix = reshape(tmp(:,2),3,[]);
zPix = reshape(tmp(:,3),3,[]);

ind1 = zeros(opts.nRowPix*opts.nColPix,1);
for ii = 1:size(sample.F,1)
ind1 = ind1|inpolygon(rPDn(2,:).',rPDn(3,:).',yPix(:,ii).',zPix(:,ii).');
end
xProjection(ind1) = 1;

end