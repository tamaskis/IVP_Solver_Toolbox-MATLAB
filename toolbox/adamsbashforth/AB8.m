%==========================================================================
%
% AB8  Propagates the state vector forward one time step using the
% Adams-Bashforth 8th-order method.
%
%   F_next = AB8(f,t,F,h)
%
% See also AB2, AB3, AB4, AB5, AB6, AB7.
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
%==========================================================================
function F = AB8(f,t,F,h)
    
    % updates F matrix with new function evaluation
    F = [F(:,2:8),f(t,F(:,9)),F(:,9)];
    
    % state vector propagated to next sample time
    F(:,9) = F(:,9)+(h/120960)*(434241*F(:,8)-1152169*F(:,7)+2183877*...
        F(:,6)-2664477*F(:,5)+2102243*F(:,4)-1041723*F(:,3)+295767*F(:,...
        2)-36799*F(:,1));
    
end