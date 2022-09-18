%==========================================================================
%
% ABM7  Propagates the state vector forward one time step using the
% Adams-Bashforth-Moulton 7th-order method.
%
%   F_next = ABM7(f,t,F,h)
%
% See also ABM2, ABM3, ABM4, ABM5, ABM6, ABM8.
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
%   F       - (p×8 double) F matrix for current sample time
%               --> columns 1-7: function evaluations at previous 7 sample
%                                times
%               -->    column 8: state (i.e. solution) at current sample 
%                                time, t
%   h       - (1×1 double) step size
%
% -------
% OUTPUT:
% -------
%   F       - (p×8 double) F matrix updated for next sample time
%               -->    column 1: function evaluation at current sample
%                                time, t
%               --> columns 2-7: function evaluations at previous 6 sample
%                                times
%               -->    column 8: state (i.e. solution) at next sample 
%                                time, t+h
%
%==========================================================================
function F = ABM7(f,t,F,h)
    
    % updates F matrix with new function evaluation
    F = [F(:,2:7),f(t,F(:,8)),F(:,8)];
    
    % predictor step
    yp = F(:,8)+(h/60480)*(198721*F(:,7)-447288*F(:,6)+705549*F(:,5)-...
        688256*F(:,4)+407139*F(:,3)-134472*F(:,2)+19087*F(:,1));
    
    % corrector step
    F(:,8) = F(:,8)+(h/60480)*(19087*f(t+h,yp)+65112*F(:,7)-46461*F(:,...
        6)+37504*F(:,5)-20211*F(:,4)+6312*F(:,3)-863*F(:,2));
    
end