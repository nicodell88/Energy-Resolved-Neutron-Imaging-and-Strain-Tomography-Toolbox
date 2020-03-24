function processFitsFiles(opts)
%PROCESSFITSFILES Processes fits files and saves 1 .mat file per projection.
%	processFitsFiles(opts)
%   Inputs:
%       - opts is a structure containing
%           opts.path		:	Directory path to fits files, can contain wildcards (*).
%           opts.filespec  	:	File name format, including type specifiers for 
%								projection number and time-of-flight index.
%           opts.proj_idx  	:	A 1-by-N vector of projection numbers to be processed.
%           opts.tof_idx    :   A 1-by-N vector of ... Only wavelengths corresponding 
%           opts.save_dir   :  	Relative or absolute file path to directory where data 
%								will be saved.
%           opts.save_str   :   File name prefix for .mat files.
%   Outputs:
%		- None.
%
% Copyright (C) 2020 The University of Newcastle, Australia
% Authors:
%   Nicholas O'Dell <Nicholas.Odell@newcastle.edu.au>
% Last modified: 23/03/2020
% This program is licensed under GNU GPLv3, see LICENSE for more details.

TBdir = fileparts(mfilename('fullpath'));
addpath(fullfile(TBdir,'utility_functions'));
addpath(fullfile(TBdir,'FitsFileProcessing'));

assert(isfield(opts,'path'),	'A file path to the data must be speciefied.')
assert(isfield(opts,'filespec'),'A file path to the data must be speciefied.')
assert(isfield(opts,'proj_idx'),'projection indicies must be specified')
assert(isfield(opts,'tof_idx'),	'projection indicies must be specified')

if ~isfield(opts,'save_dir')
	opts.save_dir = fullfile('.','preprocessed');
end

if(7~=exist(opts.save_dir,'dir'))
	mkdir(opts.save_dir)
end

delete(findall(0,'tag','TMWWaitbar'));
wh      = waitbar(0,'Processing',...
        'Name', 'Fits File Progress Bar');

for pp = 1:length(opts.proj_idx)
	
	% Get file names
	fun = @(x)sprintf(opts.filespec,opts.proj_idx(pp),x);
	fileNames = arrayfun(fun,opts.tof_idx,'UniformOutput',false);
	% Process projection
	msg = sprintf('Processing projection %d of %d',pp,length(opts.proj_idx));
	waitbar(pp/length(opts.proj_idx),wh,msg)
	data = processProjection(fileNames,opts.path,opts);
	% Save data
	msg = sprintf('Saving projection %d of %d',pp,length(opts.proj_idx));
	waitbar(pp/length(opts.proj_idx),wh,msg)
	save(fullfile(opts.save_dir,sprintf('%s_proj_%03d',opts.save_str,opts.proj_idx(pp))),'-struct','data','-v7.3')
	
end
delete(wh)

end