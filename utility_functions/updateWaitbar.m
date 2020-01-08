function [wh,flag] = updateWaitbar(n, N, wh)
%UPDATEWAITBAR Updates a waitbar
%   Initial use wh = updateWaitbar()
%   Use in for loop [wh,flag] = updateWaitbar(n, N, wh), where flag
%   indicated if the user has pressed cancel.
%
% Copyright (C) 2020 The University of Newcastle, Australia
% Authors:
%   Nicholas O'Dell <Nicholas.Odell@newcastle.edu.au>
% Last modified: 08/01/2020
% This program is licensed under GNU GPLv3, see LICENSE for more details.

flag = false;
if(nargin == 0)
    delete(findall(0,'tag','TMWWaitbar'));
    msg     = 'Fitting Bragg Edges';
    wh      = waitbar(0,msg, ...
        'Name', 'Bragg Edge Progress Bar', ...
        'CreateCancelBtn', 'setappdata(gcbf,''cancelling'',1)');
else
    complete = n/N;
    nupdates = 100;
    modu     = ceil(N/nupdates);
    
    if(mod(n, modu)==0)
        waitbar(complete,wh); % Update waitbar
    end
    
    if getappdata(wh,'cancelling') % Check if waitbar cancel button clicked
        delete(wh);
        flag = true;
    end
    
end
end

