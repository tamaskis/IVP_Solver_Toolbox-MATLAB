%==========================================================================
%
% ABM8  Propagates the state vector forward one time step using the
% Adams-Bashforth-Moulton 8th-order method.
%
%   F_next = ABM8(f,t,F,h)
%
% See also ABM2, ABM3, ABM4, ABM5, ABM6, ABM7.
%
% Copyright © 2021 Tamas Kis
% Last Update: 2022-06-07
% Website: https://tamaskis.github.io
% Contact: tamas.a.kis@outlook.com
%
% TOOLBOX DOCUMENTATION:
% https://tamaskis.github.io/IVP_Solver_Toolbox-MATLAB/
%
% TECHNICAL DOCUMENTATION:
% https://tamaskis.github.io/documentation/Solving_Initial_Value_Problems_for_ODEs.pdf
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   f       - (1×1 function_handle) dy/dt = f(t,y) --> multivariate, 
%             vector-valued function (f : ℝ×ℝᵖ → ℝᵖ) defining ODE
%   t       - (1×1 double) current sample time
%   F       - (p×9 double) F matrix for current sample time
%               --> columns 1-8: function evaluations at previous 8 sample
%                                times
%               -->    column 9: state (i.e. solution) at current sample 
%                                time, t
%   h       - (1×1 double) step size
%
% -------
% OUTPUT:
% -------
%   F       - (p×9 double) F matrix updated for next sample time
%               -->    column 1: function evaluation at current sample
%                                time, t
%               --> columns 2-8: function evaluations at previous 7 sample
%                                times
%               -->    column 9: state (i.e. solution) at next sample 
%                                time, t+h
%
% -----
% NOTE:
% -----
%   --> p = dimension of state vector (for the scalar case, p = 1)
%
%==========================================================================
function F = ABM8(f,t,F,h)
    
    % updates F matrix with new function evaluation
    F = [F(:,2:8),f(t,F(:,9)),F(:,9)];
    
    % predictor step
    yp = F(:,9)+(h/120960)*(434241*F(:,8)-1152169*F(:,7)+2183877*F(:,6)-...
        2664477*F(:,5)+2102243*F(:,4)-1041723*F(:,3)+295767*F(:,2)-...
        36799*F(:,1));
    
    % corrector step
    F(:,9) = F(:,9)+(h/120960)*(36799*f(t+h,yp)+139849*F(:,8)-121797*...
        F(:,7)+123133*F(:,6)-88547*F(:,5)+41499*F(:,4)-11351*F(:,3)+...
        1375*F(:,2));
    
end