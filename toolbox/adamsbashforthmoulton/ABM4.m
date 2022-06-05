%==========================================================================
%
% ABM4  Propagates the state vector forward one time step using the
% Adams-Bashforth-Moulton 4th-order method.
%
%   F_next = ABM4(f,t,F,h)
%
% See also ABM2, ABM3, ABM5, ABM6, ABM7, ABM8.
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
%   f       - (1×1 function_handle) dy/dt = f(t,y) --> multivariate, 
%             vector-valued function (f : ℝ×ℝᵖ → ℝᵖ) defining ODE
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
% -----
% NOTE:
% -----
%   --> p = dimension of state vector (for the scalar case, p = 1)
%
%==========================================================================
function F = ABM4(f,t,F,h)
    
    % updates F matrix with new function evaluation
    F = [F(:,2:4),f(t,F(:,5)),F(:,5)];
    
    % predictor step
    yp = F(:,5)+(h/24)*(55*F(:,4)-59*F(:,3)+37*F(:,2)-9*F(:,1));
    
    % corrector step
    F(:,5) = F(:,5)+(h/24)*(9*f(t+h,yp)+19*F(:,4)-5*F(:,3)+F(:,2));
    
end