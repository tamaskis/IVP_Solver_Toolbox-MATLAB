%==========================================================================
%
% ABM6  Propagates the state vector forward one time step using the
% Adams-Bashforth-Moulton 6th-order method.
%
%   F_next = ABM6(f,t,F,h)
%
% See also ABM2, ABM3, ABM4, ABM5, ABM7, ABM8.
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
%   F       - (p×7 double) F matrix for current sample time
%               --> columns 1-6: function evaluations at previous 6 sample
%                                times
%               -->    column 7: state (i.e. solution) at current sample 
%                                time, t
%   h       - (1×1 double) step size
%
% -------
% OUTPUT:
% -------
%   F       - (p×7 double) F matrix updated for next sample time
%               -->    column 1: function evaluation at current sample
%                                time, t
%               --> columns 2-6: function evaluations at previous 5 sample
%                                times
%               -->    column 7: state (i.e. solution) at next sample 
%                                time, t+h
%
%==========================================================================
function F = ABM6(f,t,F,h)
    
    % updates F matrix with new function evaluation
    F = [F(:,2:6),f(t,F(:,7)),F(:,7)];
    
    % predictor step
    yp = F(:,7)+(h/1440)*(4277*F(:,6)-7923*F(:,5)+9982*F(:,4)-7298*F(:,...
        3)+2877*F(:,2)-475*F(:,1));
    
    % corrector step
    F(:,7) = F(:,7)+(h/1440)*(475*f(t+h,yp)+1427*F(:,6)-798*F(:,5)+482*...
        F(:,4)-173*F(:,3)+27*F(:,2));
    
end