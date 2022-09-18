%==========================================================================
%
% ABM5  Propagates the state vector forward one time step using the
% Adams-Bashforth-Moulton 5th-order method.
%
%   F_next = ABM5(f,t,F,h)
%
% See also ABM2, ABM3, ABM4, ABM6, ABM7, ABM8.
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
%   F       - (p×6 double) F matrix for current sample time
%               --> columns 1-5: function evaluations at previous 5 sample
%                                times
%               -->    column 6: state (i.e. solution) at current sample 
%                                time, t
%   h       - (1×1 double) step size
%
% -------
% OUTPUT:
% -------
%   F       - (p×6 double) F matrix updated for next sample time
%               -->    column 1: function evaluation at current sample
%                                time, t
%               --> columns 2-5: function evaluations at previous 4 sample
%                                times
%               -->    column 6: state (i.e. solution) at next sample 
%                                time, t+h
%
%==========================================================================
function F = ABM5(f,t,F,h)
    
    % updates F matrix with new function evaluation
    F = [F(:,2:5),f(t,F(:,6)),F(:,6)];
    
    % predictor step
    yp = F(:,6)+(h/720)*(1901*F(:,5)-2774*F(:,4)+2616*F(:,3)-1274*F(:,...
        2)+251*F(:,1));
    
    % corrector step
    F(:,6) = F(:,6)+(h/720)*(251*f(t+h,yp)+646*F(:,5)-264*F(:,4)+106*...
        F(:,3)-19*F(:,2));
    
end