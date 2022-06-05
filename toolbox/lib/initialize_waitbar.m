%==========================================================================
%
% initialize_waitbar  Initializes a waitbar.
%
%   [wb,prop] = initialize_waitbar(msg)
%
% Copyright © 2022 Tamas Kis
% Last Update: 2022-06-04
% Website: https://tamaskis.github.io
% Contact: tamas.a.kis@outlook.com
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   msg     - (char) message to be displayed on waitbar
%
% -------
% OUTPUT:
% -------
%   wb      - (1×1 Figure) waitbar
%   prop    - (1×1 double) cutoff proportion to trigger waitbar update
%
%==========================================================================
function [wb,prop] = initialize_waitbar(msg)
    
    % initialize cutoff proportion needed to trigger waitbar update to 0.1
    prop = 0.1;
    
    % initializes the waitbar
    wb = waitbar(0,msg);
    
end