%==========================================================================
%
% AB4  Propagates the state vector forward one time step using the
% Adams-Bashforth 4th-order method.
%
%   F_next = AB4(f,t,F,h)
%
% See also AB2, AB3, AB5, AB6, AB7, AB8.
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
%   F       - (p×5 double) F matrix for current sample time
%               --> columns 1-4: function evaluations at previous 4 sample
%                                times
%               -->    column 5: state (i.e. solution) at current sample 
%                                time, t
%   h       - (1×1 double) step size
%
% -------
% OUTPUT:
% -------
%   F       - (p×5 double) F matrix updated for next sample time
%               -->    column 1: function evaluation at current sample
%                                time, t
%               --> columns 2-4: function evaluations at previous 3 sample
%                                times
%               -->    column 5: state (i.e. solution) at next sample 
%                                time, t+h
%
%==========================================================================
function F = AB4(f,t,F,h)
    
    % updates F matrix with new function evaluation
    F = [F(:,2:4),f(t,F(:,5)),F(:,5)];
    
    % state vector propagated to next sample time
    F(:,5) = F(:,5)+(h/24)*(55*F(:,4)-59*F(:,3)+37*F(:,2)-9*F(:,1));
    
end