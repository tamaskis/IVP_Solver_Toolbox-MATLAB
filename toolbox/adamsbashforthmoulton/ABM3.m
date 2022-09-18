%==========================================================================
%
% ABM3  Propagates the state vector forward one time step using the
% Adams-Bashforth-Moulton 3rd-order method.
%
%   F_next = ABM3(f,t,F,h)
%
% See also ABM2, ABM4, ABM5, ABM6, ABM7, ABM8.
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
%   F       - (p×4 double) F matrix for current sample time
%               --> columns 1-3: function evaluations at previous 3 sample
%                                times
%               -->    column 4: state (i.e. solution) at current sample 
%                                time, t
%   h       - (1×1 double) step size
%
% -------
% OUTPUT:
% -------
%   F       - (p×4 double) F matrix updated for next sample time
%               -->    column 1: function evaluation at current sample
%                                time, t
%               --> columns 2-3: function evaluations at previous 2 sample
%                                times
%               -->    column 4: state (i.e. solution) at next sample 
%                                time, t+h
%
%==========================================================================
function F = ABM3(f,t,F,h)
    
    % updates F matrix with new function evaluation
    F = [F(:,2:3),f(t,F(:,4)),F(:,4)];
    
    % predictor step
    yp = F(:,4)+(h/12)*(23*F(:,3)-16*F(:,2)+5*F(:,1));
    
    % corrector step
    F(:,4) = F(:,4)+(h/12)*(5*f(t+h,yp)+8*F(:,3)-F(:,2));
    
end