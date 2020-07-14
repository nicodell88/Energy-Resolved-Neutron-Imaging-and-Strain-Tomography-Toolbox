function R = stage_rotation(about_z,about_y)
% R = stage_rotation(about_z,about_y)
% builds the rotation matrix that transforms a vector from stage coordinate 
% to beamcoordinate system.
% the stage consists of two rotations, the first rotation happens about the
% y axis and the second rotation about the z axis
%
% This function is for the goniometer used in: Hendriks, J. N., Gregg, A. W.
% T., Jackson, R. R., Wensrich, C. M., Wills, A., Tremsin, A. S., ... &
% Kirstein, O. (2019). Tomographic reconstruction of triaxial strain fields
% from Bragg-edge neutron imaging. Physical Review Materials, 3(11),
% 113803.

% Copyright (C) 2020 The University of Newcastle, Australia
% Authors:
%   Nicholas O'Dell <Nicholas.Odell@newcastle.edu.au>
%   Johannes Hendriks <Johannes.Hendriks@newcastle.edu.au>
% Last modified: 20/05/2020
% This program is licensed under GNU GPLv3, see LICENSE for more details.

theta   = about_y;   % Pitch angle
psi     = about_z;   % Yaw angle

Rz = [ ...
     cos(psi), -sin(psi), 0; ...
     sin(psi),  cos(psi), 0; ...
            0,         0, 1 ...
    ];

Ry = [ ...
     cos(theta), 0, sin(theta); ...
              0, 1,          0; ...
    -sin(theta), 0, cos(theta) ...
    ];


R = Ry*Rz;

