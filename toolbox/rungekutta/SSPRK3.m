%==========================================================================
%
% SSPRK3  Propagates the state vector forward one time step using 
% the strong stability preserving Runge-Kutta third-order method.
%
%   y_next = SSPRK3(f,t,y,h)
%
% See also RK1_euler, RK2, RK2_heun, RK2_ralston, RK3, RK3_heun, 
% RK3_ralston, RK4, RK4_ralston, RK4_38.
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
%   y       - (p×1 double) state (i.e. solution) at current sample time, t
%   h       - (1×1 double) step size
%
% -------
% OUTPUT:
% -------
%   y_next  - (p×1 double) state (i.e. solution) at next sample time, t+h
%
% -----
% NOTE:
% -----
%   --> The documentation above is written specifically for the case of
%       vector-valued ODEs, but this function can also be used for matrix-
%       valued ODEs of the form dM/dt = F(t,M), where 
%       F : ℝ×ℝᵖˣʳ → ℝᵖˣʳ.
%
%==========================================================================
function y_next = SSPRK3(f,t,y,h)
    
    % k terms
    k1 = f(t,y);
    k2 = f(t+h,y+h*k1);
    k3 = f(t+h/2,y+h*k1/4+h*k2/4);
    
    % state vector propogated to next sample time
    y_next = y+(h/6)*(k1+k2+4*k3);
    
end