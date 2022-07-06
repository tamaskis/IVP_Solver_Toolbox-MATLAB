%==========================================================================
%
% update_waitbar  Updates a waitbar.
%
%   prop = update_waitbar(k,N,wb,prop)
%
% Copyright © 2022 Tamas Kis
% Last Update: 2022-07-05
% Website: https://tamaskis.github.io
% Contact: tamas.a.kis@outlook.com
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   k       - (1×1 double) current sample number (i.e. iteration)
%   N       - (1×1 double) total number of samples (i.e. iterations)
%   wb      - (1×1 Figure) waitbar
%   prop    - (1×1 double) cutoff proportion to trigger waitbar update
%
% -------
% OUTPUT:
% -------
%   prop    - (1×1 double) cutoff proportion to trigger waitbar update
%
% -----
% NOTE:
% -----
%   --> "prop" is an integer multiple of 0.1 so that the waitbar is
%       only updated after every additional 10% of progress.
%
%==========================================================================
function prop = update_waitbar(k,N,wb,prop)
    
    % only updates waitbar if current proportion exceeds cutoff proportion
    if k/N > prop
        
        % updates waitbar
        waitbar(k/N,wb);
        
        % updates cutoff proportion needed to trigger waitbar update
        prop = prop+0.1;
        
    end
    
end