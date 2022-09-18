%==========================================================================
%
% AB2  Propagates the state vector forward one time step using the
% Adams-Bashforth 2nd-order method.
%
%   F_next = AB2(f,t,F,h)
%
% See also AB3, AB4, AB5, AB6, AB7, AB8.
%
% Copyright © 2021 Tamas Kis
% Last Update: 2022-09-18
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
%   f       - (1×1 function_handle) multivariate, vector-valued function 
%             defining vector-valued ODE, dy/dt = f(t,y) (f : ℝ×ℝᵖ → ℝᵖ)
%   t       - (1×1 double) current sample time
%   F       - (p×3 double) F matrix for current sample time
%               --> columns 1-2: function evaluations at previous 2 sample
%                                times
%               -->    column 3: state (i.e. solution) at current sample 
%                                time, t
%   h       - (1×1 double) step size
%
% -------
% OUTPUT:
% -------
%   F       - (p×3 double) F matrix updated for next sample time
%               --> column 1: function evaluation at current sample time, t
%               --> column 2: function evaluation at previous sample time
%               --> column 3: state (i.e. solution) at next sample time, 
%                             t+h
%
%==========================================================================
function F = AB2(f,t,F,h)
    
    % updates F matrix with new function evaluation
    F = [F(:,2),f(t,F(:,3)),F(:,3)];
    
    % state vector propagated to next sample time
    F(:,3) = F(:,3)+(h/2)*(3*F(:,2)-F(:,1));
    
end