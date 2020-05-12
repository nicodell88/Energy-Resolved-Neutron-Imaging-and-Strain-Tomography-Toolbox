function [f,g,h] = AlignmentCostFunction(x,rPDd,beta,rVOo,rHSs,Rsh,opts)
% [f,g,h] = scanMatchCost(x,rPDd,beta,rVOo,Rns1)
% TODO

% Not solving this problem anymore
% Theta_nd    = x(1:3);     %Defines Rnd, the rotation between the detector coordinates and the beam coordinate system.
% [Rnd,dRndDX,dRndDY,dRndDZ]      = eulerRotation(Theta_nd);

% Copyright (C) 2020 The University of Newcastle, Australia
% Authors:
%   Nicholas O'Dell <Nicholas.Odell@newcastle.edu.au>
% Last modified: 12/05/2020
% This program is licensed under GNU GPLv3, see LICENSE for more details.
sigma = opts.sigma;

Theta_ns   = x(1:3);       %Defines Rns1, allows for misalignment between the goniometer and the beam.
Theta_ho   = x(4:6);       %Defines Rso1, the rotation between the sample coordinates and the stage coordinates.
rSDn        = [0;x(7:8)];   %vector to S from D in n. Stage placement.
rOHh        = x(9:11);      %vector to O from S in s. Sample placement.



[Rns,dRnsDX,dRnsDY,dRnsDZ]  = eulerRotation(Theta_ns);
[Rho,dRhoDX,dRhoDY,dRhoDZ]  = eulerRotation(Theta_ho);

nProj = length(rPDd);
nVert = length(rVOo);
nx = length(x);
%% Initialise outputs
f = zeros(nProj,1);
g = zeros(nx,nProj);
h = zeros(nx,nx,nProj);

%% Double loop
parfor i = 1:nProj
    rPDn = rPDd{i};
    
    %     Npix = length(beta{i});
    
    for j = 1:nVert
        
        %% Calculate projections onto plane
%         rVDn = (Rns*Rsh(:,:,i)* (Rho*rVOo(:,j) + rOHh) + rSDn);
        rVDn = rSDn + Rns*rHSs(:,i) + Rns*Rsh(:,:,i)*rOHh + Rns*Rsh(:,:,i)*Rho*rVOo(:,j);
        %% Calculate norm
        err = rVDn(2:3,:) - rPDn(2:3,:);
        %% Calculate cost
        ErrorSquared = dot(err,err,1);
        
        ftmp = 1/nProj * 1/(sqrt(2*pi*sigma^2)) * beta{i}.'* exp(-1/(2*sigma^2) * ErrorSquared.');
        f(i) = f(i)  + ftmp;
        
        
        
        if nargout(@AlignmentCostFunction)>1
            
            Je = [...
                dRnsDX(2:3,:)*(Rsh(:,:,i)*(Rho*rVOo(:,j) + rOHh) + rHSs(:,i)),...    dfdTheta_ns1
                dRnsDY(2:3,:)*(Rsh(:,:,i)*(Rho*rVOo(:,j) + rOHh) + rHSs(:,i)),...    dfdTheta_ns1
                dRnsDZ(2:3,:)*(Rsh(:,:,i)*(Rho*rVOo(:,j) + rOHh) + rHSs(:,i)),...    dfdTheta_ns1
                ...
                (Rns(2:3,:)* Rsh(:,:,i)*(dRhoDX*rVOo(:,j))),...  dfdTheta_so1
                (Rns(2:3,:)* Rsh(:,:,i)*(dRhoDY*rVOo(:,j))),...  dfdTheta_so1
                (Rns(2:3,:)* Rsh(:,:,i)*(dRhoDZ*rVOo(:,j))),...  dfdTheta_so1
                ...
                (eye(2)),... dfdrSDn
                ...
                (Rns(2:3,:)*Rsh(:,:,i)),... dfdrOSs
                ];
            
            
            g(:,i) = g(:,i) + ...
                1/nProj* 1/(sqrt(2*pi*sigma^2))* (-1/sigma^2) * sum(...
                Je.'*err .* beta{i}.' .* exp(-1/(2*sigma^2) * ErrorSquared)...
                ,2);
            
            %             h(:,:,i) = h(:,:,i)+...
            %                 Je.'*Je * ...
            %                 1/nProj* 1/(sqrt(2*pi*sigma^2))*  (-1/sigma^2) * ...
            %                 sum(beta{i}.' .* exp(-1/(2*sigma^2) * ErrorSquared),...
            %                 'all');
            J = (-1/sigma^2 * 1/2 * ((Je.'*err).*sqrt(beta{i}).') .* exp(-1/(4*sigma^2) * ErrorSquared)).'; 
            h(:,:,i) = h(:,:,i) + ...
                1/nProj* 2/(sqrt(2*pi*sigma^2)) * (J.'*J);
        
        end
    end
end

f = -sum(f,'all');

if nargout>1
    g =   -sum(g,2);
else
    g = [];
end

h = -sum(h,3);

% h = [];
end

