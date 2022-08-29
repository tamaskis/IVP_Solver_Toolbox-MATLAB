%==========================================================================
%
% expand_solution_arrays  Expands the arrays storing the IVP solution. 
%
%   [t_new,y_new] = expand_solution_arrays(t,y)
%
% Copyright © 2021 Tamas Kis
% Last Update: 2022-08-28
% Website: https://tamaskis.github.io
% Contact: tamas.a.kis@outlook.com
%
% TOOLBOX DOCUMENTATION:
% https://tamaskis.github.io/IVP_Solver_Toolbox-MATLAB/
%
% TECHNICAL DOCUMENTATION:
% https://tamaskis.github.io/files/Solving_Initial_Value_Problems_for_ODEs.pdf
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   t       - ((N+1)×1 double) time vector
%   y       - (p×(N+1) double) solution matrix
%
% -------
% OUTPUT:
% -------
%   t_new   - (2(N+1)×1 double) expanded time vector
%   y_new   - (p×2(N+1) double) expanded solution matrix
%
% -----
% NOTE:
% -----
%   --> This function is used to preallocate more space for the solution 
%       arrays (i.e. the time vector and the solution matrix) once they 
%       have become full.
%
%==========================================================================
function [t_new,y_new] = expand_ivp_arrays(t,y)
    
    % number of subintervals
    N = length(t)-1;
    
    % state dimension
    p = size(y,1);
    
    % expands time vector
    t_new = [t;zeros(N+1,1)];
    
    % expands solution matrix
    y_new = [y,zeros(p,N+1)];
    
end