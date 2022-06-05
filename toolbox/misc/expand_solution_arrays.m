%==========================================================================
%
% expand_solution_arrays  Expands the arrays storing the ODE solution. 
%
%   [t_new,y_new] = expand_solution_arrays(t,y)
%
% Copyright © 2021 Tamas Kis
% Last Update: 2022-06-04
% Website: https://tamaskis.github.io
% Contact: tamas.a.kis@outlook.com
%
% TOOLBOX DOCUMENTATION:
% https://tamaskis.github.io/ODE_Solver_Toolbox-MATLAB/
%
% TECHNICAL DOCUMENTATION:
% https://tamaskis.github.io/documentation/Fixed_Step_ODE_Solvers.pdf
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
%   --> p = dimension of state vector (for the scalar case, p = 1) 
%   --> N+1 = original length of time vector
%   --> This function is used to preallocate more space for the solution 
%       arrays (i.e. the time vector and the solution matrix) once they 
%       have become full.
%
%==========================================================================
function [t_new,y_new] = expand_solution_arrays(t,y)
    t_new = [t;zeros(length(t),1)];
    y_new = [y,zeros(size(y,1),length(t))];
end