%==========================================================================
%
% AB7  Propagates the state vector forward one time step using the
% Adams-Bashforth 7th-order method.
%
%   F_next = AB7(f,t,F,h)
%
% See also AB2, AB3, AB4, AB5, AB6, AB8.
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
function F = AB7(f,t,F,h)
    
    % updates F matrix with new function evaluation
    F = [F(:,2:7),f(t,F(:,8)),F(:,8)];
    
    % state vector propagated to next sample time
    F(:,8) = F(:,8)+(h/60480)*(198721*F(:,7)-447288*F(:,6)+705549*F(:,...
        5)-688256*F(:,4)+407139*F(:,3)-134472*F(:,2)+19087*F(:,1));
    
end