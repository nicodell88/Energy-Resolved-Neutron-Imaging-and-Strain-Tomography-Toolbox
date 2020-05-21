function stop = plotProj(x,state,sample,edges,rHSs,Rsh,opts,vars)
% Copyright (C) 2020 The University of Newcastle, Australia
% Authors:
%   Nicholas O'Dell <Nicholas.Odell@newcastle.edu.au>
% Last modified: 20/05/2020
% This program is licensed under GNU GPLv3, see LICENSE for more details.
stop = false;

if strcmp('iter',vars)
    
    Xopt = x;
    
    Theta_ns   = Xopt(1:3);         %Allows for misalignment between the stage and the beam
    Theta_ho    = Xopt(4:6);        %Defines Rso, the rotation between the sample coordinates and the stage coordinates.
    rSDn        = [0;Xopt(7:8)];    %vector to S from D in n. Stage placement.
    rOHh        = Xopt(9:11);       %vector to O from S in s. Sample placement.
    
    Rns = eulerRotation(Theta_ns);
    Rho = eulerRotation(Theta_ho);
    
    N = prod(opts.plotGrid);
    
    idx = sort(randsample(size(Rsh,3),N));
    
    persistent p Hdetector H
    if isempty(p) || isempty(Hdetector)
        p = gobjects(N,1);
        Hdetector = gobjects(N,1);
        H = figure(2);
        clf
        
        for i = 1:N
            subplot(opts.plotGrid(1),opts.plotGrid(2),i);
            
            %         rVDn = [rSDn] + Rns*Rs1s(:,:,idx(i))*rOSs + Rns*Rs1s(:,:,idx(i))*Rho*Ro1o(:,:,idx(i))*sample.V;
            rVDn = rSDn + Rns*rHSs(:,idx(i)) + Rns*Rsh(:,:,idx(i))*rOHh + Rns*Rsh(:,:,idx(i))*Rho*sample.V.';
            
            p(i) = patch('Faces',sample.F,'Vertices',rVDn.');
            p(i).FaceAlpha =0.2;
            axis equal
            hold on
            
            xImage = [0.02,0.02;0.02,0.02];
            zImage = -([512 512; 1 1]-513/2)*55e-6;
            yImage = -([512 1; 512 1]-513/2)*55e-6;
            
            Hdetector(i) = surf(xImage,yImage,zImage,...    %# Plot the surface
                'CData',sign(edges(:,:,idx(i))),...
                'FaceColor','texturemap');
            view([-90,0])
            msg = sprintf('Projection %d',idx(i));
            title(msg)
            %         xlim([-0.015 0.015]);
            ylim([-0.015 0.015]);
            zlim([-0.015 0.015]);
        end
    else
        
        for i = 1:N
            
            rVDn = rSDn + Rns*rHSs(:,idx(i)) + Rns*Rsh(:,:,idx(i))*rOHh + Rns*Rsh(:,:,idx(i))*Rho*sample.V.';
            p(i).Vertices = rVDn.';
            
            Hdetector(i).CData = sign(edges(:,:,idx(i)));
            msg = sprintf('Projection %d',idx(i));
            
            H.Children(N+1-i).Title.String = msg;
            %         xlim([-0.015 0.015]);
            %         ylim([-0.015 0.015]);
        end
    end
    drawnow
%     msg = sprintf('ResultsIter_%02d',state.iteration);
%     save(msg,'Theta_ns','Theta_ho','rSDn','rOHh');
%     saveas(gcf,msg,'png')
end
end