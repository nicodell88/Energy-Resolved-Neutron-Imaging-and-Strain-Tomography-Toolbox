%
% runProcessFitsFiles.m 
%
% Copyright (C) 2020 The University of Newcastle, Australia
% Authors:
%   Nicholas O'Dell <Nicholas.Odell@newcastle.edu.au>
% Last modified: 21/04/2020
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.

restoredefaultpath
clc
clear
close all

% Toolbox is designed to process fits files with the filename convention "NAME_%03d_%05d.fits"
%% Options
if isunix
options.path		= '/Volumes/Drive 1/04_Projections/*/Corrected/'; %Path including wildcards if needed
options.filespec	= 'SteelCube_%03d_%05d.fits';	%File name format including type specifiers for indexing projections and time-of-flight
options.proj_idx	= 1:2;			%Only these projections will be processed
options.tof_idx 	= 1307:1310;	%Only wavelengths corresponding to these indicies will be processed

options.save_dir	= '../preprocessed';
options.save_str    = 'steelcube';
else
options.path		= 'D:/JPARC JAN 2020/04_Projections/*/Corrected/'; %Path including wildcards if needed
options.filespec	= 'SteelCube_%03d_%05d.fits';	%File name format including type specifiers for indexing projections and time-of-flight
options.proj_idx	= 1:2;			%Only these projections will be processed
options.tof_idx 	= 1307:1310;	%Only wavelengths corresponding to these indicies will be processed

options.save_dir	= '../preprocessed';
options.save_str    = 'steelcube';
end
%% Process Data
processFitsFiles(options)



